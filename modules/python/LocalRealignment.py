from build import FRIDAY
from modules.python.ActiveRegionFinder import ActiveRegionFinder, ActiveRegionOptions
from modules.python.DeBruijnHaplotyper import DeBruijnHaplotyper


class CandidateFinderOptions(object):
    # base and map quality
    MIN_BASE_QUALITY = 15
    MIN_MAP_QUALITY = 15


class AlingerOptions(object):
    # base and map quality
    ALIGNMENT_SAFE_BASES = 20
    MIN_MAP_QUALITY = 20


class RegionBasedHaplotypes:
    def __init__(self, haplotypes, region_start, region_end):
        self.haplotypes = haplotypes
        self.region_start = region_start
        self.region_end = region_end
        self.min_read_start = None
        self.max_read_end = None
        self.reads = []

    def assign_read(self, read):
        if self.min_read_start is None or self.max_read_end is None:
            self.min_read_start = read.pos
            self.max_read_end = read.pos_end

        self.min_read_start = min(self.min_read_start, read.pos)
        self.max_read_end = max(self.max_read_end, read.pos_end)
        self.reads.append(read)

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))


class LocalAssembler:
    def __init__(self, bam_file, fasta_file, chromosome_name, region_start, region_end):
        self.bam_file_path = bam_file
        self.fasta_file_path = fasta_file
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    def perform_local_alignment(self, region_with_reads):
        if not region_with_reads.reads:
            return []
        ref_start = min(region_with_reads.min_read_start, region_with_reads.region_start) - AlingerOptions.ALIGNMENT_SAFE_BASES
        ref_end = max(region_with_reads.max_read_end, region_with_reads.region_end) + AlingerOptions.ALIGNMENT_SAFE_BASES
        fasta_handler = FRIDAY.FASTA_handler(self.fasta_file_path)

        if ref_end <= region_with_reads.region_end:
            return region_with_reads.reads
        else:
            ref_suffix = fasta_handler.get_reference_sequence(self.chromosome_name,
                                                              region_with_reads.region_end,
                                                              ref_end)

        ref_prefix = fasta_handler.get_reference_sequence(self.chromosome_name,
                                                          ref_start,
                                                          region_with_reads.region_start)
        ref = fasta_handler.get_reference_sequence(self.chromosome_name,
                                                   region_with_reads.region_start,
                                                   region_with_reads.region_end)
        ref_seq = ref_prefix + ref + ref_suffix
        haplotypes = [ref_prefix + hap + ref_suffix for hap in region_with_reads.haplotypes]

        aligner = FRIDAY.ReadAligner(ref_start, ref_end, ref_seq)

        haplotypes = sorted(set(haplotypes))

        if not haplotypes or haplotypes == [ref_seq]:
            return region_with_reads.reads

        realigned_reads = aligner.align_reads(haplotypes, region_with_reads.reads)

        return realigned_reads

    def perform_local_assembly(self, perform_alignment=True):
        # get the reads from the bam file
        bam_handler = FRIDAY.BAM_handler(self.bam_file_path)
        all_reads = bam_handler.get_reads(self.chromosome_name,
                                          self.region_start_position,
                                          self.region_end_position,
                                          CandidateFinderOptions.MIN_MAP_QUALITY,
                                          0)

        if perform_alignment is False:
            return all_reads

        # find active regions
        active_region_finder = ActiveRegionFinder(self.bam_file_path,
                                                  self.fasta_file_path,
                                                  self.chromosome_name,
                                                  self.region_start_position,
                                                  self.region_end_position)
        # mapq threshold is checked in this
        active_regions = active_region_finder.find_active_region()

        assembly_active_regions = []
        possible_regions = []

        for active_region in active_regions:
            start_pos, end_pos = active_region

            if end_pos - start_pos > ActiveRegionOptions.MAX_REGION_SIZE:
                continue

            db_graph = DeBruijnHaplotyper(self.bam_file_path,
                                          self.fasta_file_path,
                                          self.chromosome_name,
                                          start_pos,
                                          end_pos)

            reference_sequence, haplotypes = db_graph.find_haplotypes()

            if haplotypes:
                assembly_active_regions.append(RegionBasedHaplotypes(haplotypes, start_pos, end_pos))
                possible_regions.append((start_pos, end_pos))

        if not possible_regions:
            return all_reads
        # now we have the list that we filtered at the beginning of this script
        realigned_reads = list()
        for read in all_reads:
            read_start = read.pos
            read_end = read.pos_end
            read_range = (read_start, read_end)
            overlapping_lengths = [RegionBasedHaplotypes.overlap_length_between_ranges(region, read_range)
                                   for region in possible_regions]
            max_length = max(overlapping_lengths)
            if max_length <= 0:
                realigned_reads.append(read)
                continue

            max_window_index = max(range(len(possible_regions)), key=lambda i: overlapping_lengths[i])
            assembly_active_regions[max_window_index].assign_read(read)

        for active_region in assembly_active_regions:
            realigned_reads.extend(self.perform_local_alignment(active_region))

        return realigned_reads
