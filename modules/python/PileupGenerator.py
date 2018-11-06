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

    def decode_image_row(self, read_array):
        global_base_color_reverse = {250: 'A', 30: 'C', 180: 'G', 100: 'T', 0: ' ', 10: 'N'}
        for pixel in read_array:
            base = pixel[0]
            print(global_base_color_reverse[base], end='')
        print()

    def generate_pileup(self, reads, windows, positional_candidates):
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
        pileup_images = image_generator.create_window_pileups(windows, reads)
        # print("GOT READS")

        # for pileup_image in pileup_images:
        #     print(pileup_image.chromosome_name, pileup_image.start_pos, pileup_image.end_pos)
        #     for image_row in pileup_image.image:
        #         self.decode_image_row(image_row)
        #
        # exit()
        return pileup_images
        # exit(0)
        # for read in reads:
        #     if read.read_id not in read_window_map:
        #         continue
        #     img_row, read_pos = image_generator.read_to_image_row(read)
            # print(read_pos, img_row)

