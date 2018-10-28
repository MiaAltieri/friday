import argparse
import sys
import operator
import os
import torch
from torch.utils.data import DataLoader
from torchvision import transforms
from modules.python.dataloader_fake import PileupDataset, TextColor
from collections import defaultdict
from modules.python.vcf_writer import VCFWriter
from tqdm import tqdm


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


def predict(test_file, batch_size, num_workers):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param test_file: File to predict on
    :param batch_size: Batch size used for prediction
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # the prediction table/dictionary
    prediction_dict = defaultdict(list)
    transformations = transforms.Compose([transforms.ToTensor()])

    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    # data loader
    validation_data = PileupDataset(test_file, transformations)
    validation_loader = DataLoader(validation_data,
                                   batch_size=batch_size,
                                   shuffle=False,
                                   num_workers=num_workers
                                   )

    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    with torch.no_grad():
        for labels, records in tqdm(validation_loader, file=sys.stdout, ncols=50):
            # record each of the predictions from a batch prediction
            for i in range(0, labels.size(0)):
                rec = records[i]
                chr_name, pos_st, pos_end, ref, alt1, alt2, rec_type_alt1, rec_type_alt2 = rec.rstrip().split('\t')[0:8]

                true_label = labels[i]
                fake_probs = [0.0] * 3
                fake_probs[true_label] = 1.0

                prob_hom, prob_het, prob_hom_alt = fake_probs
                prediction_dict[pos_st].append((chr_name, pos_st, pos_end, ref, alt1, alt2, rec_type_alt1, rec_type_alt2,
                                                prob_hom, prob_het, prob_hom_alt))

    return prediction_dict


def produce_vcf(prediction_dictionary, bam_file_path, sample_name, output_dir):
    """
    Convert prediction dictionary to a VCF file
    :param prediction_dictionary: prediction dictionary containing predictions of each image records
    :param bam_file_path: Path to the BAM file
    :param sample_name: Name of the sample in the BAM file
    :param output_dir: Output directory
    :return:
    """
    # object that can write and handle VCF
    vcf_writer = VCFWriter(bam_file_path, sample_name, output_dir)

    # print("test5")

    # collate multi-allelic records to a single record
    all_calls = []
    for pos in sorted(prediction_dictionary.keys()):
        records = prediction_dictionary[pos]
        # if record is multi-allelic then combine and get the genotype
        if len(records) > 1:
            chrm, st_pos, end_pos, ref, alts, genotype, qual, gq = vcf_writer.get_genotype_for_multiple_allele(records)
        else:
            chrm, st_pos, end_pos, ref, alts, genotype, qual, gq = vcf_writer.get_genotype_for_single_allele(records)

        all_calls.append((chrm, int(st_pos), int(end_pos), ref, alts, genotype, qual, gq))

    # sort based on position
    all_calls.sort(key=operator.itemgetter(1))
    last_end = 0
    for record in all_calls:
        # get the record filter ('PASS' or not)
        rec_filter = vcf_writer.get_filter(record, last_end)
        # get proper alleles. INDEL alleles are handled here.
        record = vcf_writer.get_proper_alleles(record)
        chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual, phred_gq = record
        # if genotype is not HOM keep track of where the previous record ended
        if genotype != '0/0':
            last_end = end_pos
        # add the record to VCF
        vcf_writer.write_vcf_record(chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual, phred_gq, rec_filter)


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

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--test_file",
        type=str,
        required=True,
        help="Testing data description csv file."
    )
    parser.add_argument(
        "--bam_file",
        type=str,
        required=True,
        help="Path to the BAM file."
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
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "SAMPLE NAME: " + FLAGS.sample_name + "\n")
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PLEASE USE --sample_name TO CHANGE SAMPLE NAME.\n")
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "OUTPUT DIRECTORY: " + FLAGS.output_dir + "\n")
    record_dict = predict(FLAGS.test_file, FLAGS.batch_size, FLAGS.num_workers)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "PREDICTION COMPLETED SUCCESSFULLY.\n")
    produce_vcf(record_dict, FLAGS.bam_file, FLAGS.sample_name, FLAGS.output_dir)
    sys.stderr.write(TextColor.GREEN + "INFO: " + TextColor.END + "FINISHED CALLING VARIANT.\n")
    sys.stderr.flush()
