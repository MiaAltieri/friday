"""
            possible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom       hom     no
            het       hom     yes
            het       het     yes
            hom_alt   hom     yes
            hom       None    no
            het       None    yes
            hom_alt   None    yes

            impossible combos:
            gt1       gt2     Candidate validated?
            --------------------------------------
            hom_alt   hom_alt NA
            het       hom_alt NA
            hom       hom_alt NA
            hom_alt   het     NA
            hom       het     NA
            None      None    NA

"""
GENOTYPE_DICT = {"Hom": 0, "Het": 1, "Hom_alt": 2}
PLOIDY = 2

# DEBUG_FREQUENCIES = False
DEBUG_PRINT_ALL = False

# Candidate data indexes
CHR_NAME = 0
START = 1
STOP = 2
REF_INDEX = 3
ALT1 = 4
ALT2 = 5
ALT1_TYPE = 6
ALT2_TYPE = 7

# Positional vcf indexes
SNP_CANDIDATE, IN_CANDIDATE, DEL_CANDIDATE = 1, 2, 3
SNP_DEL = 3

# VCF record indexes
REF, ALT, GT = 0, 1, 2

# Genotype codes
HOM, HET, HOM_ALT = 0, 1, 2
from build import FRIDAY


class CandidateLabeler:
    def __init__(self, vcf_path, chromosome_name, start_pos, end_pos):
        """
        Initialize candidateLabeler object
        :param fasta_handler: module that fetches reference sequence substrings from a FASTA file
        """
        self.vcf_offset = -1                    # pysam vcf coords are 1-based ... >:[ ... this is what Kishwar wanted
        self.delete_char = '*'
        vcf_handler = FRIDAY.VCF_handler(vcf_path)
        self.positional_vcf = vcf_handler.get_positional_vcf_records(chromosome_name, start_pos, end_pos)

    def _handle_insert(self, rec):
        """
        Process a record that has an insert
        :param rec: VCF record
        :return: attributes of the record
        """
        ref_seq = rec.ref   # no change necessary
        alt_seq = rec.alt   # no change necessary

        pos = rec.pos + self.vcf_offset
        return pos, ref_seq, alt_seq, rec.type

    def _handle_delete(self, rec):
        """
        Process a record that has deletes.
        Deletes are usually grouped together, so we break each of the deletes to make a list.
        :param rec: VCF record containing a delete
        :return: A list of delete attributes
        """
        delete_list = []
        for i in range(0, len(rec.ref)):
            if i < len(rec.alt):
                continue
            ref_seq = rec.ref[i]
            alt_seq = '*'
            pos = rec.pos + i + self.vcf_offset
            genotype = rec.type
            delete_list.append((pos, ref_seq, alt_seq, genotype))
        return delete_list

    @staticmethod
    def _resolve_suffix_for_insert(ref, alt):
        len_ref = len(ref) - 1
        if len_ref == 0:
            return ref, alt
        suffix_removed_alt = alt[:-len_ref]
        return ref[0], suffix_removed_alt

    @staticmethod
    def _resolve_suffix_for_delete(ref, alt):
        len_alt = len(alt) - 1
        if len_alt == 0:
            return ref, alt
        suffix_removed_ref = ref[:-len_alt]
        return suffix_removed_ref, alt[0]

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
    def _is_supported(genotypes):
        """
        Check if genotype has anything other than Hom
        :param genotypes: Genotype tuple
        :return: Boolean [True if it has Het of Hom_alt]
        """
        supported = False
        gt_set = set(genotypes)

        if len(gt_set) == 0:
            supported = False
        elif len(gt_set) == 1:
            if HOM in gt_set:
                supported = False
            else:
                supported = True
        elif len(gt_set) > 1:
            supported = True

        return supported

    def _is_position_supported(self, genotypes):
        """
        Check if a position has any genotype other than Hom
        :param genotypes: Genotypes list of that position
        :return: Boolean [True if it has Het or Hom_alt]
        """
        in_supported = self._is_supported(genotypes[IN])
        del_supported = self._is_supported(genotypes[DEL])
        snp_supported = self._is_supported(genotypes[SNP])

        return in_supported or del_supported or snp_supported

    def _generate_list(self, chromosome_name, start, stop, alleles_snp, alleles_in, ref_seq, genotypes):
        """
        Generate a list of attributes that can be saved of a labeled candidate
        :param chromosome_name: Name of chromosome
        :param start: Allele start position
        :param stop: Allele end position
        :param alleles_snp: SNP alleles
        :param alleles_in: Insert alleles
        :param ref_seq: reference Sequence
        :param genotypes: Genotypes
        :return: A list containing (chr start stop ref_seq alt1 alt2 gt1 gt2)
        """
        all_candidates = []
        for i, allele_tuple in enumerate(alleles_snp):
            allele, freq = allele_tuple
            gt = genotypes[SNP][i][0]
            gt_q = genotypes[SNP][i][1]
            gt_f = genotypes[SNP][i][2]
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, gt, gt_q, gt_f])

        for i, allele_tuple in enumerate(alleles_in):
            allele, freq = allele_tuple
            gt = genotypes[IN][i][0]
            gt_q = genotypes[IN][i][1]
            gt_f = genotypes[IN][i][2]
            all_candidates.append([chromosome_name, start, stop, ref_seq, allele, gt, gt_q, gt_f])

        return all_candidates

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
            return 0, 0
        if gt1 == HET and gt2 == HOM:
            return 1, 0
        if gt1 == HOM and gt2 == HET:
            return 2, 0
        if gt1 == HET and gt2 == HET:
            return 1, 2
        if gt1 == HOM_ALT and gt2 == HOM:
            return 1, 1
        if gt2 == HOM_ALT and gt2 == HOM:
            return 2, 2
        else:
            import sys
            sys.stderr.write("ERROR: GENOTYPES DIDNOT MATCH: ", gt1, gt2)
        return None

    def get_labeled_candidates(self, candidate_sites):
        """
        Label candidates given variants from a VCF
        :param positional_vcf: IN/DEL/SNPs separated into VCF records, expanded into a 1-to-1 ref_pos:variant allele
        :param candidate_sites: Candidates
        :return: List of labeled candidate sites
        """
        # list of all labeled candidates
        all_labeled_candidates = []
        # for candidate in candidate_sites:
        #     candidate.print()
        # print("-----------------")
        # print(type(self.positional_vcf))
        # for pos in list(self.positional_vcf.keys()):
        #     print(self.positional_vcf[pos])
        #     for rec in self.positional_vcf[pos]:
        #         print(rec.chromosome_name, rec.start_pos, end=' ')
        #         for alt in rec.alt_allele:
        #             print(alt.ref, alt.alt_allele, alt.alt_type)
        #     print(self.positional_vcf[pos])
        # print("--------------")
        # for each candidate
        for candidate_site in candidate_sites:
            # test the alleles across IN, DEL, and SNP variant dictionaries
            genotypes = self.get_combined_gt(self.get_label_of_allele(candidate=candidate_site))
            candidate_site.set_genotype(genotypes)
            all_labeled_candidates.append(candidate_site)
            # candidate_site.print()
            # exit()
            # # get a list of attributes that can be saved
            # labeled_candidate = candidate_site + genotypes
            # all_labeled_candidates.append(labeled_candidate)
            #
            # if DEBUG_PRINT_ALL:
            #     print()
            #     print(allele_start, allele_stop, ref_sequence)
            #     print("alleles:        ", alleles)
            #     print("Genotypes:     ", genotypes)

        return all_labeled_candidates
