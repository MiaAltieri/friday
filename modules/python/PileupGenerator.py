from build import FRIDAY
from collections import defaultdict
from modules.python.Options import CandidateFinderOptions


class PileupGenerator:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
        self.chromosome_name = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def generate_pileup(self, reads, windows, positional_candidates):
        # find out which read goes to what window. Reads that do not overlap with any window will not be processed
        read_window_map = defaultdict(list)
        for read_id in range(len(reads)):
            reads[read_id].set_read_id(read_id)
            for window_id in range(len(windows)):
                if self.overlap_length_between_ranges((reads[read_id].pos, reads[read_id].pos_end),
                                                      windows[window_id]):
                    read_window_map[read_id].append(window_id)

        ref_start = max(0, windows[0][0] - CandidateFinderOptions.SAFE_BASES)
        ref_end = windows[-1][1] + CandidateFinderOptions.SAFE_BASES
        reference_sequence = self.fasta_handler.get_reference_sequence(self.chromosome_name,
                                                                       ref_start,
                                                                       ref_end)
        # image generator object
        image_generator = FRIDAY.ImageGenerator(reference_sequence,
                                                self.chromosome_name,
                                                ref_start,
                                                ref_end,
                                                positional_candidates)

        for read in reads:
            if read.read_id not in read_window_map:
                continue
            img_row, read_pos = image_generator.read_to_image_row(read)
            # print(read_pos, img_row)
