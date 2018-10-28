from build import FRIDAY


class CandidateFinderOptions(object):
    MIN_MAPPING_QUALITY = 14
    MIN_BASE_QUALITY = 15
    SAFE_BASES = 50
    # the linear regression model is used inside C++ code


class CandidateFinder:
    def __init__(self, fasta_file_path, contig, start, end):
        self.fasta_file_path = fasta_file_path
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def find_candidates(self, reads):
        # ref_start = self.region_start
        # ref_end = self.region_end
        # for read in reads:
        #     ref_start = min(ref_start, read.pos)
        #     ref_end = max(ref_end, read.pos_end)
        # get the reference from the fasta file
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)
        reference_sequence = fasta_handler.get_reference_sequence(self.contig,
                                                                  self.region_start - CandidateFinderOptions.SAFE_BASES,
                                                                  self.region_end + CandidateFinderOptions.SAFE_BASES)
        # find the active region
        candidate_finder = FRIDAY.CandidateFinder(reference_sequence,
                                                  self.contig,
                                                  self.region_start,
                                                  self.region_end,
                                                  self.region_start - CandidateFinderOptions.SAFE_BASES,
                                                  self.region_end + CandidateFinderOptions.SAFE_BASES)

        # find active regions
        candidates = candidate_finder.find_candidates(reads)

        # print("RAW CANDIDATES")
        # for candidate in candidates:
        #     candidate.print()
        # print("--------------------")
        # exit()
        # candidates_in_region = []
        # for candidate in candidates:
        #     if(self.overlap_length_between_ranges((candidate.pos, candidate.pos_end),
        #                                           (self.region_start, self.region_end))):
        #         candidates_in_region.append(candidate)

        return candidates
