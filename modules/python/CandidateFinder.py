from build import FRIDAY


class CandidateFinderOptions(object):
    MIN_MAPPING_QUALITY = 15
    MIN_BASE_QUALITY = 15
    # the linear regression model is used inside C++ code


class CandidateFinder:
    def __init__(self, fasta_file_path, contig, start, end):
        self.fasta_file_path = fasta_file_path
        self.contig = contig
        self.region_start = start
        self.region_end = end

    def find_candidates(self, reads):
        # get the reference from the fasta file
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)
        reference_sequence = fasta_handler.get_reference_sequence(self.contig, self.region_start, self.region_end)
        # find the active region
        candidate_finder = FRIDAY.CandidateFinder(reference_sequence,
                                                  self.contig,
                                                  self.region_start,
                                                  self.region_end)

        # find active regions
        candidates = candidate_finder.find_candidates(reads)

        for candidate in candidates:
            print(candidate.chromosome_name, candidate.pos, candidate.pos_end, candidate.ref, candidate.alt1, candidate.alt2, candidate.alt1_type, candidate.alt2_type)
        return candidates
