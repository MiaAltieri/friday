import argparse
import math
import time
import os
import sys
import pickle
import h5py
import numpy as np

from build import FRIDAY
from modules.python.IntervalTree import IntervalTree
from modules.python.LocalRealignment import LocalAssembler
from modules.python.CandidateFinder import CandidateFinder
from modules.python.TextColor import TextColor
from modules.python.TsvHandler import TsvHandler
from modules.python.FileManager import FileManager
from modules.python.PileupGenerator import PileupGenerator
from modules.python.Options import ImageSizeOptions
from modules.python.CandidateLabler import CandidateLabeler
from modules.python.AlignmentSummarizer import AlignmentSummarizer
from modules.python.DataStore import DataStore
"""
This script creates training images from BAM, Reference FASTA and truth VCF file. The process is:
- Find candidates that can be variants
- Label candidates using the VCF
- Create images for each candidate

Input:
- BAM file: Alignment of a genome
- REF file: The reference FASTA file used in the alignment
- VCF file: A truth VCF file
- BED file: A confident bed file. If confident_bed is passed it will only generate train set for those region.

Output:
- TENSOR files: Containing images and their labels.
- CSV file: Containing records of images and their location in the tensor file.
"""

# Global debug helpers
DEBUG_PRINT_CANDIDATES = False
DEBUG_TIME_PROFILE = False
DEBUG_TEST_PARALLEL = False
BED_POSITION_BUFFER = 0


class View:
    """
    Process manager that runs sequence of processes to generate images and their labebls.
    """
    def __init__(self, chromosome_name, bam_file_path, reference_file_path, vcf_path, train_mode, confident_tree, downsample_rate):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param reference_file_path: Path to the reference FASTA file
        :param vcf_path: Path to the VCF file
        :param confident_tree: Dictionary containing all confident trees. NULL if parameter not passed.
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = reference_file_path
        self.vcf_path = vcf_path
        self.bam_handler = FRIDAY.BAM_handler(bam_file_path)
        self.fasta_handler = FRIDAY.FASTA_handler(reference_file_path)
        self.train_mode = train_mode
        self.confident_tree = confident_tree[chromosome_name] if confident_tree else None
        self.interval_tree = IntervalTree(self.confident_tree) if confident_tree else None
        self.downsample_rate = downsample_rate

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    @staticmethod
    def build_chromosomal_interval_trees(confident_bed_path):
        """
        Produce a dictionary of intervals trees, with one tree per chromosome
        :param confident_bed_path: Path to confident bed file
        :return: trees_chromosomal
        """
        # create an object for tsv file handling
        tsv_handler_reference = TsvHandler(tsv_file_path=confident_bed_path)
        # create intervals based on chromosome
        intervals_chromosomal_reference = tsv_handler_reference.get_bed_intervals_by_chromosome(universal_offset=-1)

        return intervals_chromosomal_reference

    def get_labeled_candidate_sites(self, selected_candidate_list, start_pos, end_pos, filter_hom_ref=False):
        """
        Lable selected candidates of a region and return a list of records
        :param selected_candidate_list: List of all selected candidates with their alleles
        :param start_pos: start position of the region
        :param end_pos: end position of the region
        :param filter_hom_ref: whether to ignore hom_ref VCF records during candidate validation
        :return: labeled_sites: Labeled candidate sites. Each containing proper genotype.
        """
        candidate_labler = CandidateLabeler(self.fasta_path, self.vcf_path)
        labeled_candidates = candidate_labler.get_labeled_candidates(self.chromosome_name,
                                                                     start_pos,
                                                                     end_pos,
                                                                     selected_candidate_list)
        del candidate_labler
        return labeled_candidates

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    @staticmethod
    def a_fully_contains_range_b(range_a, range_b):
        if range_b[0] >= range_a[0] and range_b[1] <= range_a[1]: return True
        return False

    def parse_region(self, start_position, end_position, local_alignment_flag):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
        """
        # st_time = time.time()
        # print("STARTING: ", self.chromosome_name, start_position, end_position)
        confident_intervals_in_region = []
        if self.train_mode:
            confident_intervals_in_region = self.interval_tree.find(start_position, end_position)
            if not confident_intervals_in_region:
                return 0, 0, None, None

        local_assembler = LocalAssembler(self.bam_handler,
                                         self.fasta_handler,
                                         self.chromosome_name,
                                         start_position,
                                         end_position)

        reads = local_assembler.perform_local_assembly(self.downsample_rate, perform_alignment=local_alignment_flag)

        # returning here to see if local assembly is happening or not
        if not reads:
            return 0, 0, None, None

        candidate_finder = CandidateFinder(self.fasta_handler,
                                           self.chromosome_name,
                                           start_position,
                                           end_position)
        candidate_list, position_to_read_map = candidate_finder.find_candidates(reads)

        if not candidate_list:
            return len(reads), 0, None, None

        image_generator = PileupGenerator(self.fasta_handler,
                                          self.chromosome_name,
                                          start_position,
                                          end_position)

        # # get all labeled candidate sites
        if self.train_mode:
            if not candidate_list:
                return 0, 0, None, None

            # should summarize and get labels here
            labeled_candidates = self.get_labeled_candidate_sites(candidate_list, start_position, end_position,
                                                                  True)

            confident_candidates = []
            # for subsetting the candidates
            for candidate in labeled_candidates:
                for interval in confident_intervals_in_region:
                    if self.a_fully_contains_range_b(interval, (candidate.pos_start, candidate.pos_start)):
                        confident_candidates.append(candidate)

            labeled_candidates = confident_candidates

            candidates,  images = image_generator.generate_pileup(reads, labeled_candidates, position_to_read_map,
                                                                  train_mode=True)

            return len(reads), len(candidates), candidates, images
        else:
            if not candidate_list:
                return 0, 0, None, None
            candidates, images = image_generator.generate_pileup(reads, candidate_list, position_to_read_map,
                                                                 train_mode=False)
            return len(reads), len(candidates), candidates, images


def create_output_dir_for_chromosome(output_dir, chr_name):
    """
    Create an internal directory inside the output directory to dump choromosomal summary files
    :param output_dir: Path to output directory
    :param chr_name: chromosome name
    :return: New directory path
    """
    path_to_dir = output_dir + chr_name + "/"
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    summary_path = path_to_dir + "summary" + "/"
    if not os.path.exists(summary_path):
        os.mkdir(summary_path)

    return path_to_dir


def chromosome_level_parallelization(chr_list,
                                     bam_file,
                                     ref_file,
                                     vcf_file,
                                     confident_intervals,
                                     output_path,
                                     total_threads,
                                     thread_id,
                                     train_mode,
                                     downsample_rate,
                                     local_alignment,
                                     max_size=1000):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_list: List of chromosomes to be processed
    :param bam_file: path to BAM file
    :param ref_file: path to reference FASTA file
    :param vcf_file: path to VCF file
    :param max_size: Maximum size of a segment
    :param output_path: path to output directory
    :return:
    """
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = FRIDAY.FASTA_handler(ref_file)
    data_file_name = output_path + "images" + "_thread_" + str(thread_id) + ".hdf"

    data_file = DataStore(data_file_name, mode='w')

    for chr_name, region in chr_list:
        if not region:
            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) + 1)
        else:
            interval_start, interval_end = tuple(region)

        all_intervals = []
        for pos in range(interval_start, interval_end, max_size):
            all_intervals.append((pos, min(interval_end, pos + max_size - 1)))

        intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

        view = View(chromosome_name=chr_name,
                    bam_file_path=bam_file,
                    reference_file_path=ref_file,
                    vcf_path=vcf_file,
                    train_mode=train_mode,
                    confident_tree=confident_intervals,
                    downsample_rate=downsample_rate)

        start_time = time.time()
        total_reads_processed = 0
        total_candidates = 0

        for count, interval in enumerate(intervals):
            _start, _end = interval
            n_reads, n_candidates, candidates, images = view.parse_region(start_position=_start,
                                                                          end_position=_end,
                                                                          local_alignment_flag=local_alignment)

            total_reads_processed += n_reads
            total_candidates += n_candidates

            if not candidates:
                continue

            candidate_data_file_name = output_path + "candidates" + "_thread_" + str(thread_id) + \
                                       "_" + str(interval[0]) + "_" + str(interval[1]) + ".pkl"

            with open(candidate_data_file_name, 'wb') as f:
                pickle.dump(candidates, f, pickle.HIGHEST_PROTOCOL)


            # data_file.write_images(images, chr_name)
            # data_file.write_candidates(candidates, chr_name)

        print("CHROMOSOME: ", chr_name,
              "THREAD ID: ", thread_id,
              "READS: ", total_reads_processed,
              "CANDIDATES: ", total_candidates,
              "TOTAL TIME ELAPSED: ", int(math.floor(time.time()-start_time)/60), "MINS",
              time.time()-start_time, "SEC")


def summary_file_to_csv(output_dir_path, chr_list):
    """
    Remove the abundant number of summary files and bind them to one
    :param output_dir_path: Path to the output directory
    :param chr_list: List of chromosomes
    :return:
    """
    for chr_name in chr_list:
        # here we dumped all the bed files
        path_to_dir = output_dir_path + chr_name + "/summary/"

        concatenated_file_name = output_dir_path + chr_name + ".csv"

        filemanager_object = FileManager()
        # get all bed file paths from the directory
        file_paths = filemanager_object.get_file_paths_from_directory(path_to_dir)
        # dump all bed files into one
        filemanager_object.concatenate_files(file_paths, concatenated_file_name)
        # delete all temporary files
        filemanager_object.delete_files(file_paths)
        # remove the directory
        os.rmdir(path_to_dir)


def handle_output_directory(output_dir, thread_id):
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


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s:
    :return:
    """
    if s.lower() not in {'false', 'true'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true'


def get_chromosme_list(chromosome_names):
    split_names = chromosome_names.strip().split(',')
    split_names = [name.strip() for name in split_names]

    chromosome_name_list = []
    for name in split_names:
        # split on region
        region = None
        if ':' in name:
            name_region = name.strip().split(':')

            if len(name_region) != 2:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

            name, region = tuple(name_region)
            region = region.strip().split('-')
            region = [int(pos) for pos in region]

            if len(region) != 2 or not region[0] <= region[1]:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

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
                chromosome_name_list.append((chr_prefix + str(chr_seq), region))
        else:
            chromosome_name_list.append((name, region))

    return chromosome_name_list


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads of interest."
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        default=None,
        help="VCF file path."
    )
    parser.add_argument(
        "--bed",
        type=str,
        default=None,
        help="Path to confident BED file"
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number E.g.: 3"
    )
    parser.add_argument(
        "--train_mode",
        type=boolean_string,
        default=False,
        help="If true then a dry test is run."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--thread_id",
        type=int,
        required=False,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--downsample_rate",
        type=float,
        required=False,
        default=1.0,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--local_alignment",
        type=boolean_string,
        default=True,
        help="If true then perform local alignment."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = get_chromosme_list(FLAGS.chromosome_name)

    # if the confident bed is not empty then create the tree
    if FLAGS.bed:
        confident_intervals = View.build_chromosomal_interval_trees(FLAGS.bed)
    else:
        confident_intervals = None

    if FLAGS.train_mode and (not confident_intervals or not FLAGS.vcf):
        sys.stderr.write(TextColor.RED + "ERROR: TRAIN MODE REQUIRES --vcf AND --bed TO BE SET.\n" + TextColor.END)
        exit(1)
    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir), FLAGS.thread_id)

    chromosome_level_parallelization(chr_list,
                                     FLAGS.bam,
                                     FLAGS.fasta,
                                     FLAGS.vcf,
                                     confident_intervals,
                                     output_dir,
                                     FLAGS.threads,
                                     FLAGS.thread_id,
                                     FLAGS.train_mode,
                                     FLAGS.downsample_rate,
                                     FLAGS.local_alignment)

