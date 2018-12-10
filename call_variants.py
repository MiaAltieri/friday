import argparse
import sys
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader
from torchvision import transforms
import multiprocessing
from modules.python.models.ModelHander import ModelHandler
from modules.python.models.dataloader_predict import SequenceDataset
from modules.python.Options import ImageSizeOptions
from modules.python.TextColor import TextColor
from collections import defaultdict
from modules.python.vcf_writer import VCFWriter
from modules.python.FileManager import FileManager
import operator
import pickle
from tqdm import tqdm
import os
"""
This script uses a trained model to call variants on a given set of images generated from the genome.
The process is:
- Create a prediction table/dictionary using a trained neural network
- Convert those predictions to a VCF file

INPUT:
- A trained model
- Set of images for prediction

Output:
- A VCF file containing all the variants.
"""
SNP = 1
IN = 2
DEL = 3
HOM = 0
HET = 1
HOM_ALT = 2

prediction_dict = defaultdict(list)
reference_dict = defaultdict(tuple)


def predict(test_file, batch_size, model_path, gpu_mode, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param test_file: File to predict on
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param gpu_mode: If true, predictions will be done over GPU
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # the prediction table/dictionary
    prediction_dict.clear()
    reference_dict.clear()
    chromosome_name = ''
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    test_dset = SequenceDataset(test_file, transformations)
    testloader = DataLoader(test_dset,
                            batch_size=batch_size,
                            shuffle=False,
                            num_workers=num_workers
                            )

    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    # load the model
    encoder_model, decoder_model, hidden_size, gru_layers, prev_ite = \
        ModelHandler.load_model_for_training(model_path,
                                             input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                             seq_len=ImageSizeOptions.SEQ_LENGTH,
                                             num_classes=ImageSizeOptions.TOTAL_LABELS)
    hidden_size = 512
    gru_layers = 3

    if gpu_mode:
        encoder_model = encoder_model.cuda()
        encoder_model = torch.nn.DataParallel(encoder_model).cuda()
        decoder_model = decoder_model.cuda()
        decoder_model = torch.nn.DataParallel(decoder_model).cuda()

    # Change model to 'eval' mode (BN uses moving mean/var).
    encoder_model.eval()
    decoder_model.eval()

    sys.stderr.write(TextColor.PURPLE + 'MODEL LOADED\n' + TextColor.END)
    with torch.no_grad():
        for images, dictionary_path, positional_info in tqdm(testloader, ncols=50):
            start_index = ImageSizeOptions.CONTEXT_SIZE
            end_index = start_index + (ImageSizeOptions.SEQ_LENGTH - 2 * ImageSizeOptions.CONTEXT_SIZE)

            if gpu_mode:
                # encoder_hidden = encoder_hidden.cuda()
                images = images.cuda()

            unrolling_genomic_position = np.zeros((images.size(0)), dtype=np.int64)
            encoder_hidden = torch.FloatTensor(images.size(0), gru_layers * 2, hidden_size).zero_()
            context_vector, hidden_encoder = encoder_model(images, encoder_hidden)

            for seq_index in range(start_index, end_index):
                batches = images.size(0)

                current_batch_size = images.size(0)
                attention_index = torch.from_numpy(np.asarray([seq_index] * current_batch_size)).view(-1, 1)

                attention_index_onehot = torch.FloatTensor(current_batch_size, images.size(2))

                attention_index_onehot.zero_()
                attention_index_onehot.scatter_(1, attention_index, 1)

                output_dec, decoder_hidden, attn = decoder_model(attention_index_onehot,
                                                                 context_vector=context_vector,
                                                                 encoder_hidden=hidden_encoder)

                m = nn.Softmax(dim=1)
                soft_probs = m(output_dec)
                output_preds = soft_probs.cpu()

                for batch in range(batches):

                    allele_dict_path = dictionary_path[batch]
                    chr_name, start_pos, end_pos = positional_info[batch].split(' ')
                    # needs to be a check here
                    chromosome_name = chr_name
                    # current_genomic_position = int(start_positions[batch])
                    current_genomic_position = int(start_pos) + unrolling_genomic_position[batch]

                    preds = output_preds[batch, :].data
                    top_n, top_i = preds.topk(1)
                    predicted_label = top_i[0].item()
                    reference_dict[current_genomic_position] = allele_dict_path
                    prediction_dict[current_genomic_position].append((predicted_label, preds))

                    unrolling_genomic_position[batch] += 1

    return chromosome_name


def get_record_from_prediction(pos, positional_record):
    predictions = prediction_dict[pos]
    genotype, qual, gq = VCFWriter.process_prediction(pos, predictions)
    return positional_record, genotype, qual, gq


def produce_vcf_records(chromosome_name, output_dir, thread_no, pos_list):
    """
    Convert prediction dictionary to a VCF file
    :param: arg_tuple: Tuple of arguments containing these values:
    - chromosome_name: Chromosome name
    - pos_list: List of positions where we will search for variants
    - prediction_dict: prediction dictionary containing predictions of each image records
    - reference_dict: Dictionary containing reference information
    - bam_file_path: Path to the BAM file
    - sample_name: Name of the sample in the BAM file
    - output_dir: Output directory
    - thread_id: Unique id assigned to each thread
    :return:
    """
    # object that can write and handle VCF
    # vcf_writer = VCFWriter(bam_file_path, sample_name, output_dir, thread_id)
    # collate multi-allelic records to a single record
    current_allele_dict = ''
    allele_dict = {}
    record_file = open(output_dir + chromosome_name + "_" + str(thread_no) + ".tsv", 'w')
    for pos in pos_list:
        allele_dict_path = reference_dict[pos]
        if allele_dict_path != current_allele_dict:
            allele_dict = pickle.load(open(allele_dict_path, 'rb'))
            current_allele_dict = allele_dict_path

        if pos not in allele_dict:
            continue
        positional_record = allele_dict[pos] if pos in allele_dict else None
        if positional_record is None:
            continue

        positional_record, genotype, qual, gq = get_record_from_prediction(pos, positional_record)

        if genotype == '0/0':
            continue
        # print('BEFORE', record)
        ref, alts, genotype = VCFWriter.get_proper_alleles(positional_record, genotype)
        # print('AFTER', record)
        if len(alts) == 1:
            alts.append('.')
        rec_end = int(pos + len(ref) - 1)
        record_string = chromosome_name + "\t" + str(pos) + "\t" + str(rec_end) + "\t" + ref + "\t" + '\t'.join(alts) \
                        + "\t" + genotype + "\t" + str(qual) + "\t" + str(gq) + "\t" + "\n"
        record_file.write(record_string)


def merge_call_files(vcf_file_directory):
    filemanager_object = FileManager()
    # get all bed file paths from the directory
    file_paths = filemanager_object.get_file_paths_from_directory(vcf_file_directory)
    all_records = []
    for file_path in file_paths:
        with open(file_path, 'r') as tsv:
            for line in tsv:
                chr_name, pos_st, pos_end, ref, alt1, alt2, genotype, qual, gq = line.strip().split('\t')
                alts = []
                pos_st, pos_end, qual, gq = int(pos_st), int(pos_end), float(qual), float(gq)
                if alt1 != '.':
                    alts.append(alt1)
                if alt2 != '.':
                    alts.append(alt2)
                all_records.append((chr_name, pos_st, pos_end, ref, alts, genotype, qual, gq))

    filemanager_object.delete_files(file_paths)

    return all_records


def call_variant(csv_file, batch_size, model_path, gpu_mode, num_workers, vcf_dir, vcf_writer, max_threads):
    chr_name = predict(csv_file, batch_size, model_path, gpu_mode, num_workers)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION GENERATED SUCCESSFULLY.\n")

    pos_list = list(prediction_dict.keys())

    if not pos_list:
        sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "VARIANT CALLING COMPLETE.\n")
        return

    each_chunk_size = max(1, int(len(pos_list) / max_threads))
    thread_no = 1
    # produce_vcf_records(chr_name, vcf_dir, thread_no, pos_list)
    # exit()

    for i in tqdm(range(0, len(pos_list), each_chunk_size), file=sys.stdout, ncols=50):
        start_position = i
        end_position = min(i + each_chunk_size, len(pos_list) + 1)

        sub_pos = pos_list[start_position:end_position]
        # gather all parameters
        args = (chr_name, vcf_dir, thread_no, sub_pos)
        p = multiprocessing.Process(target=produce_vcf_records, args=args)
        p.start()
        thread_no += 1

        # wait until we have room for new processes to start
        while True:
            if len(multiprocessing.active_children()) < max_threads:
                break

    # wait until we have room for new processes to start
    while True:
        if len(multiprocessing.active_children()) == 0:
            break
    all_calls = merge_call_files(vcf_dir)

    # sort based on position
    all_calls.sort(key=operator.itemgetter(1))
    # print(all_calls)
    last_end = 0

    for record in all_calls:
        # get the record filter ('PASS' or not)
        rec_filter = VCFWriter.get_filter(record, last_end)
        # get proper alleles. INDEL alleles are handled here.
        # record = VCFWriter.get_proper_alleles(record)
        chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual, phred_gq = record
        # if genotype is not HOM keep track of where the previous record ended
        if genotype != '0/0':
            # HOM
            last_end = end_pos
        # add the record to VCF
        vcf_writer.write_vcf_record(chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual, phred_gq, rec_filter)

    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "VARIANT CALLING COMPLETE.\n")


def call_variants_on_multiple_chromosome(csv_dir, chr_list, bam_file_path, sample_name, output_dir,
                                         batch_size, model_path, gpu_mode, num_workers, max_threads):
    # object that can write and handle VCF
    vcf_writer = VCFWriter(bam_file_path, sample_name, output_dir)

    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "SAMPLE NAME: " + sample_name + "\n")
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PLEASE USE --sample_name TO CHANGE SAMPLE NAME.\n")
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + output_dir + "\n")

    for chromosome_name in chr_list:
        test_file = csv_dir + chromosome_name + ".csv"
        sys.stderr.write(TextColor.BLUE + "INFO: PREDICTING " + TextColor.END + str(chromosome_name) + "\n")
        call_variant(test_file,
                     batch_size,
                     model_path,
                     gpu_mode,
                     num_workers,
                     vcf_dir,
                     vcf_writer,
                     max_threads)
        sys.stderr.write(TextColor.BLUE + "INFO: " + TextColor.END + "PREDICTION COMPLETED SUCCESSFULLY.\n")
        sys.stderr.flush()
    os.rmdir(vcf_dir)


def get_chromosme_list(chromosome_names):
    split_names = chromosome_names.strip().split(',')
    split_names = [name.strip() for name in split_names]

    chromosome_name_list = []
    for name in split_names:
        range_split = name.split('-')
        if len(range_split) > 1:
            chr_prefix = ''
            for p in name:
                if p.isdigit():
                    break
                else:
                    chr_prefix = chr_prefix + p

            int_ranges = []
            for item in range_split:
                s = ''.join(i for i in item if i.isdigit())
                int_ranges.append(int(s))
            int_ranges = sorted(int_ranges)

            for chr_seq in range(int_ranges[0], int_ranges[-1] + 1):
                chromosome_name_list.append(chr_prefix + str(chr_seq))
        else:
            chromosome_name_list.append(name)

    return chromosome_name_list


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    vcf_path = output_dir + "vcfs" + "/"
    if not os.path.exists(vcf_path):
        os.mkdir(vcf_path)

    return output_dir, vcf_path


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv_dir",
        type=str,
        required=True,
        help="CSV directory containing all image segments for prediction."
    )
    parser.add_argument(
        "--bam_file",
        type=str,
        required=True,
        help="Path to the BAM file."
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number E.g.: 3"
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--model_path",
        type=str,
        default='./CNN.pkl',
        help="Saved model path."
    )
    parser.add_argument(
        "--gpu_mode",
        type=bool,
        default=False,
        help="If true then cuda is on."
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        required=False,
        default='NA12878',
        help="Sample name of the sequence."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default='vcf_output',
        help="Output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir, vcf_dir = handle_output_directory(FLAGS.output_dir)
    chr_list = get_chromosme_list(FLAGS.chromosome_name)

    call_variants_on_multiple_chromosome(FLAGS.csv_dir,
                                         chr_list,
                                         FLAGS.bam_file,
                                         FLAGS.sample_name,
                                         FLAGS.output_dir,
                                         FLAGS.batch_size,
                                         FLAGS.model_path,
                                         FLAGS.gpu_mode,
                                         FLAGS.num_workers,
                                         FLAGS.num_workers)
