from collections import defaultdict
import numpy as np
import torch
from build import FRIDAY
"""
This script creates pileup images given a vcf record, bam alignment file and reference fasta file.

imageChannels: Handles how many channels to create for each base and their structure

"""

SNP_CANDIDATE, IN_CANDIDATE, DEL_CANDIDATE = 1, 2, 3
# Genotype codes
HOM, HET, HOM_ALT = 0, 1, 2

REF_BAND_SIZE = 5
MAX_COLOR_VALUE = 254
BASE_QUALITY_CAP = 40
MAP_QUALITY_CAP = 60
MAP_QUALITY_FILTER = 5
MIN_DELETE_QUALITY = 0
IMAGE_DEPTH_THRESHOLD = 300

# '*'->insert '.'->delete
global_base_color = {'A': 250, 'C': 30, 'G': 180, 'T': 100, '.': 0, '*': 0, 'N': 10}
global_base_color_reverse = {250: 'A', 30: 'C', 180: 'G', 100: 'T', 0: ' ', 10: 'N'}



class CIGAR_OPERATIONS(object):
    MATCH = 0
    IN = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PAD = 6
    EQUAL = 7
    DIFF = 8
    BACK = 9
    UNSPECIFIED = -1


class ImageGenerator:
    """
    Processes a pileup around a position
    """
    def __init__(self, vcf_path, ref_seq, chromosome_name, ref_start, ref_end, candidate_map):
        self.ref_seq = ref_seq
        self.chromosome_name = chromosome_name
        self.ref_start = ref_start
        self.ref_end = ref_end

        vcf_handler = FRIDAY.VCF_handler(vcf_path)
        self.positional_vcf = vcf_handler.get_positional_vcf_records(chromosome_name, ref_start - 20, ref_end + 20)

        self.read_to_tensor = defaultdict()
        self.candidate_map = candidate_map

    def get_ref_sequence(self, start, end):
        _st = start - self.ref_start
        _end = end - self.ref_end

        return self.ref_seq[_st:_end]

    def get_which_allele(self, pos, ref, alt, alt_type):
        if pos not in self.candidate_map.keys():
            return 0

        candidate = self.candidate_map[pos]
        # candidate.print()
        if ref == candidate.ref and \
           alt == candidate.alt1 and \
           alt_type == candidate.alt1_type:
            return 1

        if ref == candidate.ref and \
           alt == candidate.alt2 and \
           alt_type == candidate.alt2_type:
            return 2

        return 0

    def convert_read_to_tensor(self, read):
        cigartuples = read.cigar_tuples
        read_sequence = read.sequence
        ref_pos = read.pos
        read_index = 0
        map_qual_color = int(MAX_COLOR_VALUE * (min(read.mapping_quality, MAP_QUALITY_CAP) / MAP_QUALITY_CAP))
        strand_color = 240 if read.flags.is_reverse else 70
        read_array = []

        _st = max(min(read.pos, self.ref_end), self.ref_start)
        # _end = min(read.pos_end, self.ref_end)

        for i, cigar in enumerate(cigartuples):
            cigar_op = cigar.cigar_op
            cigar_len = cigar.cigar_len
            if cigar_op == CIGAR_OPERATIONS.MATCH or \
               cigar_op == CIGAR_OPERATIONS.DIFF or \
               cigar_op == CIGAR_OPERATIONS.EQUAL:
                cigar_index = 0
                if ref_pos < self.ref_start:
                    cigar_index = min(self.ref_start - ref_pos, cigar_len)
                    read_index += cigar_index
                    ref_pos += cigar_index

                if ref_pos >= self.ref_end:
                    break

                for i in range(cigar_index, cigar_len):
                    read_allele = read_sequence[read_index]
                    ref_allele = self.ref_seq[ref_pos - self.ref_start]

                    alt_channel = 0
                    if read_allele != ref_allele:
                        alt_channel = self.get_which_allele(ref_pos, ref_allele, read_allele, SNP_CANDIDATE)
                        # print("SNP", read.read_id, ref_pos, read_allele, ref_allele, alt_channel)

                    base_quality_color = int(MAX_COLOR_VALUE * (min(read.base_qualities[read_index], BASE_QUALITY_CAP) / BASE_QUALITY_CAP))
                    base_color = global_base_color[read_allele] if read_allele in global_base_color else 0
                    alt_color = 0
                    if alt_channel != 0:
                        alt_color = 240 if alt_channel == 1 else 125

                    read_pixel = [base_color, base_quality_color, map_qual_color, strand_color, alt_color]
                    read_array.append(read_pixel)

                    ref_pos += 1
                    read_index += 1
                    if ref_pos >= self.ref_end:
                        break

            if cigar_op == CIGAR_OPERATIONS.IN:
                if self.ref_start <= ref_pos - 1 < self.ref_end:
                    if read_array:
                        anchor_pixel = read_array.pop()
                    if i == 0:
                        _st -= 1
                    read_allele = read_sequence[read_index:read_index+cigar_len]
                    ref_allele = self.ref_seq[ref_pos-self.ref_start - 1]
                    in_allele = ref_allele + read_allele
                    alt_channel = self.get_which_allele(ref_pos - 1, ref_allele, in_allele, IN_CANDIDATE)
                    # print("IN: ", ref_pos - 1, ref_allele, in_allele, alt_channel)

                    # pixel construction
                    base_color = global_base_color['*'] if '*' in global_base_color else 0
                    avg_base_quality = sum(read.base_qualities[read_index:read_index+cigar_len]) / cigar_len

                    base_quality_color = int(MAX_COLOR_VALUE * (min(avg_base_quality, BASE_QUALITY_CAP) / BASE_QUALITY_CAP))

                    alt_color = 0
                    if alt_channel != 0:
                        alt_color = 240 if alt_channel == 1 else 125

                    anchor_pixel = [base_color, base_quality_color, map_qual_color, strand_color, alt_color]
                    read_array.append(anchor_pixel)

                read_index += cigar_len

            if cigar_op == CIGAR_OPERATIONS.DEL or \
               cigar_op == CIGAR_OPERATIONS.REF_SKIP or \
               cigar_op == CIGAR_OPERATIONS.PAD:
                if self.ref_start <= ref_pos - 1 < self.ref_end:
                    ref_allele = self.ref_seq[ref_pos-self.ref_start - 1]
                    read_allele = self.ref_seq[ref_pos-self.ref_start - 1:ref_pos-self.ref_start + cigar_len]
                    alt_channel = self.get_which_allele(ref_pos - 1, ref_allele, read_allele, DEL_CANDIDATE)
                    # print("DEL: ", ref_pos - 1, ref_allele, read_allele, alt_channel)

                    if read_array:
                        anchor_pixel = read_array.pop()
                    if i == 0:
                        _st -= 1

                    # pixel construction, most of the channels will be 0 as it's delete
                    alt_color = 0
                    if alt_channel != 0:
                        alt_color = 240 if alt_channel == 1 else 125

                    anchor_pixel = [[0, 0, 0, strand_color, alt_color]] * (cigar_len + 1)
                    read_array.extend(anchor_pixel)
                    # print("DEL: ", ref_pos - 1, ref_allele, read_allele, alt_channel)
                ref_pos += cigar_len

            if cigar_op == CIGAR_OPERATIONS.SOFT_CLIP:
                # if self.ref_start <= ref_pos - 1 < self.ref_end:
                #     read_allele = read_sequence[read_index:read_index+cigar_len]
                #     ref_allele = self.ref_seq[ref_pos-self.ref_start - 1]
                #     in_allele = ref_allele + read_allele
                #     alt_channel = self.get_which_allele(ref_pos - 1, ref_allele, in_allele, IN_CANDIDATE)
                #     print("SOFT CLIP: ", ref_pos - 1, ref_allele, in_allele, alt_channel)
                read_index += cigar_len

            # if cigar_op == CIGAR_OPERATIONS.REF_SKIP or cigar_op == CIGAR_OPERATIONS.PAD:
            #     ref_pos += cigar_len

        return read_array, _st, min(ref_pos, self.ref_end)

    def decode_image_row(self, read_array):
        for pixel in read_array:
            base = pixel[0]
            print(global_base_color_reverse[base], end='')
        print()

    def get_ref_row(self, ref_seq):
        ref_array = []
        for base in ref_seq:
            base_color = global_base_color[base] if base in global_base_color else 0

            pixel = [base_color, MAX_COLOR_VALUE, MAX_COLOR_VALUE, MAX_COLOR_VALUE, MAX_COLOR_VALUE]
            ref_array.append(pixel)
        return ref_array

    def generate_image(self, window, reads, image_height=100):
        window_start, window_end = window

        ref_seq = self.ref_seq[window_start-self.ref_start:window_end-self.ref_start]
        ref_array = self.get_ref_row(ref_seq)
        # print(window_end-window_start)
        whole_image = []
        whole_image.extend([ref_array] * REF_BAND_SIZE)

        for read in reads:
            if read.pos > window_end or read.pos_end <= window_start:
                continue

            if read.read_id not in self.read_to_tensor:
                # process the read from self.read_to_tensor[read.read_id]
                self.read_to_tensor[read.read_id] = self.convert_read_to_tensor(read)

            read_array, _st, _end = self.read_to_tensor[read.read_id]
            read_segment_array = []
            segment_read_start = 0
            empties_on_left = 0
            segment_read_end = len(read_array)
            empties_on_right = 0

            if _end <= window_start:
                continue
            if _st >= window_end:
                continue

            if _st < window_start:
                segment_read_start = window_start - _st
            elif _st > window_start:
                empties_on_left = _st - window_start

            if _end > window_end:
                segment_read_end = window_end - _st
            elif _end < window_end:
                empties_on_right = window_end - _end
            left_empty = [[0, 0, 0, 0, 0]] * empties_on_left
            right_empty = [[0, 0, 0, 0, 0]] * empties_on_right
            # print(_st, _end, window_start, window_end)
            # print(segment_read_start, segment_read_end)
            if segment_read_start == segment_read_end:
                segment_read_end += 1
            core_values = read_array[segment_read_start:segment_read_end]
            read_segment_array.extend(left_empty)
            read_segment_array.extend(core_values)
            read_segment_array.extend(right_empty)

            if len(read_segment_array) != 20:
                print('READ', _st, _end, len(read_array))
                print('WINDOW', window_start, window_end)
                print(read.pos_end <= window_start)
                print(segment_read_start, segment_read_end)
                print(ref_seq)
                self.decode_image_row(read_segment_array)
                print('Window', window_start, window_end)
                print('Read', _st, _end, len(read_array), len(read_segment_array))
                print(len(left_empty), len(right_empty), len(core_values))
                print(read.pos, read.sequence)
                for ct in read.cigar_tuples:
                    print(ct.cigar_op, ct.cigar_len)
                print(segment_read_start, segment_read_end)
                print(len(left_empty), len(core_values), len(right_empty))
                exit()

            if len(whole_image) < image_height:
                whole_image.append(read_segment_array)

            if len(whole_image) == image_height:
                break

        # for row in whole_image:
        #     self.decode_image_row(row)
        # print("-------------------------")

        empty_rows = image_height - len(whole_image)
        for i in range(empty_rows):
            empty_row = [[0, 0, 0, 0, 0]] * (window_end - window_start)
            whole_image.append(empty_row)

        np_array_image = np.array(whole_image, dtype=np.uint8)
        np_array_image = np_array_image.transpose(2, 0, 1)



        return torch.from_numpy(np_array_image)

    def get_label_of_allele(self, candidate):
        """
        Given positional VCFs (IN, DEL, SNP), variant type and a candidate allele, return the try genotype.
        :param positional_vcf: Three dictionaries for each position
        :param candidate_allele: Candidate allele
        :param allele_types: Alt allele type: IN,DEL,SNP
        :return: genotype
        """
        gts = [0, 0]
        if candidate.pos not in self.positional_vcf:
            return gts

        for rec in self.positional_vcf[candidate.pos]:
            for alt in rec.alt_allele:
                # alt_ref = alt.ref
                # alt_alt = alt.alt_allele
                # if alt.alt_type == IN_CANDIDATE:
                #     alt_ref, alt_alt = self._resolve_suffix_for_insert(alt.ref, alt.alt_allele)
                # if alt.alt_type == DEL_CANDIDATE:
                #     alt_ref, alt_alt = self._resolve_suffix_for_delete(alt.ref, alt.alt_allele)

                # print(alt_ref, candidate.ref)
                # print(alt_alt, candidate.alt1)
                # print(alt.alt_type, candidate.alt1_type)
                # print(".....")
                if alt.ref == candidate.ref and \
                        alt.alt_allele == candidate.alt1 and \
                        alt.alt_type == candidate.alt1_type:
                    if rec.genotype[0] == 1 and rec.genotype[1] == 1:
                        gts[0] = HOM_ALT
                    elif rec.genotype[0] == 1 or rec.genotype[1] == 1:
                        gts[0] = HET
                # print(alt_ref, candidate.ref)
                # print(alt_alt, candidate.alt2)
                # print(alt.alt_type, candidate.alt2_type)
                if alt.ref == candidate.ref and \
                        alt.alt_allele == candidate.alt2 and \
                        alt.alt_type == candidate.alt2_type:
                    if rec.genotype[0] == 2 and rec.genotype[1] == 2:
                        gts[1] = HOM_ALT
                    elif rec.genotype[0] == 2 or rec.genotype[1] == 2:
                        gts[1] = HET
                # print("-------")
        return gts

    @staticmethod
    def get_combined_gt(gt):
        """
        Given two genotypes get the combined genotype. This is used to create labels for the third image.
        If two alleles have two different genotypes then the third genotype is inferred using this method.

        - If genotype1 is HOM then genotype of third image is genotype2
        - If genotype2 is HOM then genotype of third image is genotype1
        - If both gt are  HOM then genotype of third image is HOM
        - If genotype1, genotype2 both are HET then genotype of third image is HOM_ALT
        - If none of these cases match then we have an invalid genotype
        :param gt1: Genotype of first allele
        :param gt2: Genotype of second allele
        :return: genotype of image where both alleles are used together
        """
        gt1, gt2 = gt
        if gt1 == HOM and gt2 == HOM:
            return 0 # 0/0
        if gt1 == HET and gt2 == HOM:
            return 1 # 1/0
        if gt1 == HOM_ALT and gt2 == HOM:
            return 2 # 1/1
        if gt1 == HOM and gt2 == HET:
            return 3 # 0/2
        if gt2 == HOM_ALT and gt2 == HOM:
            return 4 # 2/2
        if gt1 == HET and gt2 == HET:
            return 5 # 1/2
        else:
            import sys
            sys.stderr.write("ERROR: GENOTYPES DIDNOT MATCH: ", gt1, gt2)
        return None

    def get_window_label(self, window):
        window_label = []
        for pos in range(window[0], window[1]):
            if pos in self.candidate_map:
                gt = self.get_combined_gt(self.get_label_of_allele(self.candidate_map[pos]))
                window_label.append(gt)
            else:
                window_label.append(0)

        return torch.from_numpy(np.array(window_label, dtype=np.uint8))

    def generate_labeled_images(self, windows, reads):
        _id = 1
        for read in reads:
            read.set_read_id(_id)
            _id += 1

        # need three things here -> window, image, label, candidate dictionary name
        generated_image_info = []
        for window in windows:
            label_tensor = self.get_window_label(window)
            image_tensor = self.generate_image(window, reads)
            generated_image_info.append(((self.chromosome_name, window[0], window[1]), image_tensor, label_tensor))

        return generated_image_info
