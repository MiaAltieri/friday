import argparse
import math
import time
import os
import sys
import multiprocessing
from tqdm import tqdm

from build import FRIDAY
from modules.python.bam_handler import BamHandler
from modules.python.LocalRealignment import LocalAssembler
from modules.python.CandidateFinder import CandidateFinder
from modules.python.CandidateLabler import CandidateLabeler
from modules.python.TextColor import TextColor
from modules.python.TsvHandler import TsvHandler
from modules.python.FileManager import FileManager
from modules.python.ImageGenerator import ImageGenerator
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
- H5PY files: Containing images and their label of the genome.
- CSV file: Containing records of images and their location in the H5PY file.
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
    def __init__(self, chromosome_name, bam_file_path, reference_file_path, vcf_path, train_mode, summary_writer):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param reference_file_path: Path to the reference FASTA file
        :param vcf_path: Path to the VCF file
        :param output_file_path: Path to the output directory where images are saved
        :param confident_tree: Dictionary containing all confident trees. NULL if parameter not passed.
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = reference_file_path
        self.vcf_path = vcf_path
        self.summary_writer = summary_writer
        self.bam_handler = FRIDAY.BAM_handler(bam_file_path)
        self.fasta_handler = FRIDAY.FASTA_handler(reference_file_path)
        self.train_mode = train_mode

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

        return labeled_candidates

    def parse_region(self, start_position, end_position):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :param thread_no: Thread no for this region
        :return:
        """
        # st_time = time.time()
        # print("STARTING", thread_no, start_position, end_position)
        local_assembler = LocalAssembler(self.bam_handler,
                                         self.fasta_handler,
                                         self.chromosome_name,
                                         start_position,
                                         end_position)

        reads = local_assembler.perform_local_assembly()

        candidate_finder = CandidateFinder(self.fasta_handler,
                                           self.chromosome_name,
                                           start_position,
                                           end_position)
        candidates = candidate_finder.find_candidates(reads)

        # # get all labeled candidate sites
        if self.train_mode:
            labeled_sites = self.get_labeled_candidate_sites(candidates, start_position, end_position, True)
        else:
            labeled_sites = candidates
        #
        # # if DEBUG_PRINT_CANDIDATES:
        # #     for candidate in labeled_sites:
        # #         print(candidate)
        #
        # # generate and save candidate images
        ImageGenerator.generate_and_save_candidate_images(labeled_sites, self.summary_writer)

        return len(reads), len(candidates)
        #
        # # end_time = time.time()
        # print("ELAPSED ", thread_no, start_position, end_position, end_time - st_time)


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


def chromosome_level_parallelization(chr_name,
                                     bam_file,
                                     ref_file,
                                     vcf_file,
                                     confident_intervals,
                                     output_path,
                                     total_threads,
                                     thread_id,
                                     train_mode,
                                     max_size=1000):
    """
    This method takes one chromosome name as parameter and chunks that chromosome in max_threads.
    :param chr_name: Name of the chromosome
    :param bam_file: path to BAM file
    :param ref_file: path to reference FASTA file
    :param vcf_file: path to VCF file
    :param max_size: Maximum size of a segment
    :param output_path: path to output directory
    :return:
    """
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = FRIDAY.FASTA_handler(ref_file)

    interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name))
    # interval_start, interval_end = (265759, 269859)

    all_intervals = []
    for pos in range(interval_start, interval_end, max_size):
        all_intervals.append((pos, pos + max_size))

    intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

    smry = None
    if intervals:
        smry = open(output_path + "summary" + '_' + chr_name + "_" + str(thread_id) + ".csv", 'w')

    view = View(chromosome_name=chr_name,
                bam_file_path=bam_file,
                reference_file_path=ref_file,
                vcf_path=vcf_file,
                summary_writer=smry,
                train_mode=train_mode)

    start_time = time.time()
    total_reads_processed = 0
    total_candidates = 0
    for interval in intervals:
        _start, _end = interval
        n_reads, n_candidates = view.parse_region(start_position=_start, end_position=_end)
        total_reads_processed += n_reads
        total_candidates += n_candidates

    print("TOTAL TIME ELAPSED: ", thread_id, time.time()-start_time,
          'READS: ', total_reads_processed,
          'CANDIDATES: ', total_candidates)


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


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s:
    :return:
    """
    if s.lower() not in {'false', 'true'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true'


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
        required=True,
        help="VCF file path."
    )
    parser.add_argument(
        "--bed",
        type=str,
        default='',
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
    FLAGS, unparsed = parser.parse_known_args()
    # if the confident bed is not empty then create the tree
    # if FLAGS.bed != '':
    #     confident_intervals = View.build_chromosomal_interval_trees(FLAGS.bed)
    # else:
    #     confident_intervals = None
    #
    # if confident_intervals is not None:
    #     sys.stderr.write(TextColor.PURPLE + "CONFIDENT TREE LOADED\n" + TextColor.END)
    # else:
    #     sys.stderr.write(TextColor.RED + "CONFIDENT BED IS NULL\n" + TextColor.END)
    confident_intervals = None
    chromosome_level_parallelization(FLAGS.chromosome_name,
                                     FLAGS.bam,
                                     FLAGS.fasta,
                                     FLAGS.vcf,
                                     confident_intervals,
                                     FLAGS.output_dir,
                                     FLAGS.threads,
                                     FLAGS.thread_id,
                                     FLAGS.train_mode)

