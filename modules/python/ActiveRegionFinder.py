from build import FRIDAY


class ActiveRegionOptions(object):
    MIN_REGION_SIZE = 80
    MAX_REGION_SIZE = 1000
    REGION_EXPANSION = 25
    # the linear regression model is used inside C++ code


class ActiveRegionFinder:
    def __init__(self, bam_file_path, fasta_file_path, contig, start, end):
        self.bam_file_path = bam_file_path
        self.fasta_file_path = fasta_file_path
        self.contig = contig
        self.region_start = max(0, start - ActiveRegionOptions.REGION_EXPANSION)
        self.region_end = end + ActiveRegionOptions.REGION_EXPANSION

    def find_active_region(self):
        # get the reads from the bam file
        bam_handler = FRIDAY.BAM_handler(self.bam_file_path)
        reads = bam_handler.get_reads(self.contig, self.region_start, self.region_end)

        # get the reference from the fasta file
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)
        reference_sequence = fasta_handler.get_reference_sequence(self.contig, self.region_start, self.region_end)

        # find the active region
        active_region_finder = FRIDAY.ActiveRegionFinder(reference_sequence,
                                                         self.contig,
                                                         self.region_start,
                                                         self.region_end)

        # find active regions
        active_regions = active_region_finder.find_active_region(reads)

        return active_regions
