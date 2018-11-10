from pysam import VariantFile, VariantHeader
from modules.python.bam_handler import BamHandler
import math
import time
import numpy as np

DEL_TYPE = 3
IN_TYPE = 2
SNP_TYPE = 1
HOM = 0
HET = 1
HOM_ALT = 2


class VCFWriter:
    def __init__(self, bam_file_path, sample_name, output_dir):
        self.bam_handler = BamHandler(bam_file_path)
        bam_file_name = bam_file_path.rstrip().split('/')[-1].split('.')[0]
        vcf_header = self.get_vcf_header(sample_name)
        time_str = time.strftime("%m%d%Y_%H%M%S")

        self.vcf_file = VariantFile(output_dir + bam_file_name + '_' + time_str + '.vcf', 'w', header=vcf_header)

    def write_vcf_record(self, chrm, st_pos, end_pos, ref, alts, genotype, qual, gq, rec_filter):
        alleles = tuple([ref]) + tuple(alts)
        genotype = self.get_genotype_tuple(genotype)
        end_pos = int(end_pos) + 1
        st_pos = int(st_pos)

        vcf_record = self.vcf_file.new_record(contig=str(chrm), start=st_pos, stop=end_pos, id='.', qual=qual,
                                              filter=rec_filter, alleles=alleles, GT=genotype, GQ=gq)
        self.vcf_file.write(vcf_record)

    @staticmethod
    def prediction_label_to_allele(label):
        label_to_allele = {0:  ['0', '0'],  1:  ['0', '1'], 2:  ['1', '1'], 3:  ['0', '2'], 4:  ['2', '2'],
                           5:  ['1', '2']}
        return label_to_allele[label]

    @staticmethod
    def get_qual_and_gq(probabilities, predicted_class):
        qual = 1.0 - probabilities[0]
        phred_qual = min(60, -10 * np.log10(1 - qual) if 1 - qual >= 0.0000001 else 60)
        phred_qual = math.ceil(phred_qual * 100.0) / 100.0

        gq = probabilities[predicted_class]
        phred_gq = min(60, -10 * np.log10(1 - gq) if 1 - gq >= 0.0000001 else 60)
        phred_gq = math.ceil(phred_gq * 100.0) / 100.0
        return phred_qual, phred_gq

    @staticmethod
    def solve_multiple_alts(alts, ref):
        type1, type2 = alts[0][1], alts[1][1]
        alt1, alt2 = alts[0][0], alts[1][0]
        if type1 == DEL_TYPE and type2 == DEL_TYPE:
            if len(alt2) > len(alt1):
                return alt2, ref, alt2[0] + alt2[len(alt1):]
            else:
                return alt1, ref, alt1[0] + alt1[len(alt2):]
        elif type1 == IN_TYPE and type2 == IN_TYPE:
            return ref, alt1, alt2
        elif type1 == DEL_TYPE or type2 == DEL_TYPE:
            if type1 == DEL_TYPE and type2 == IN_TYPE:
                return alt1, ref, alt2 + alt1[1:]
            elif type1 == IN_TYPE and type2 == DEL_TYPE:
                return alt2, alt1 + alt2[1:], ref
            elif type1 == DEL_TYPE and type2 == SNP_TYPE:
                return alt1, ref, alt2 + alt1[1:]
            elif type1 == SNP_TYPE and type2 == DEL_TYPE:
                return alt2, alt1 + alt2[1:], ref
            elif type1 == DEL_TYPE:
                return alt1, ref, alt2
            elif type2 == DEL_TYPE:
                return alt2, alt1, ref
        else:
            return ref, alt1, alt2

    @staticmethod
    def solve_single_alt(alts, ref):
        alt1, alt_type = alts[0]
        if alt_type == DEL_TYPE:
            return alt1, ref, '.'

        return ref, alt1, '.'

    @staticmethod
    def get_genotype_tuple(genotype):
        split_values = genotype.split('/')
        split_values = [int(x) for x in split_values]
        return tuple(split_values)

    @staticmethod
    def process_prediction(pos, prediction_alt1, prediction_alt2):
        # get the list of prediction labels
        # assume both are homozygous first
        alt1_probability = [0.0, 0.0, 0.0]
        alt2_probability = [0.0, 0.0, 0.0]
        if prediction_alt1:
            count = 0
            for label, probability in prediction_alt1:
                count += 1
                for j, prob_value in enumerate(probability):
                    alt1_probability[j] += prob_value
            alt1_probability = [prob / count for prob in alt1_probability]
        if prediction_alt2:
            count = 0
            for label, probability in prediction_alt2:
                count += 1
                for j, prob_value in enumerate(probability):
                    alt2_probability[j] += prob_value
            alt1_probability = [prob / count for prob in alt1_probability]
        # probability that the site genotype is 0/0
        p00 = min(alt1_probability[0], alt2_probability[0])
        p01 = alt1_probability[1]
        p11 = alt1_probability[2]
        p02 = alt2_probability[1]
        p22 = alt2_probability[2]
        p12 = min(max(alt1_probability[1], alt1_probability[2]),
                  max(alt2_probability[1], alt2_probability[2]))

        # print(alt_probs)
        prob_list = [p00, p01, p11, p02, p22, p12]
        # print(prob_list)
        sum_probs = sum(prob_list)
        # print(sum_probs)
        normalized_list = [(float(i) / sum_probs) if sum_probs else 0 for i in prob_list]
        prob_list = normalized_list
        # print(prob_list)
        # print(sum(prob_list))
        gq, index = 0, 0
        for i, prob in enumerate(prob_list):
            if gq <= prob and prob > 0:
                index = i
                gq = prob
        # get alts from label
        genotype = VCFWriter.prediction_label_to_allele(index)
        genotype = genotype[0] + '/' + genotype[1]

        qual = sum(prob_list) - prob_list[0]
        phred_qual = min(60, -10 * np.log10(1 - qual) if 1 - qual >= 0.0000001 else 60)
        phred_qual = math.ceil(phred_qual * 100.0) / 100.0
        phred_gq = min(60, -10 * np.log10(1 - gq) if 1 - gq >= 0.0000001 else 60)
        phred_gq = math.ceil(phred_gq * 100.0) / 100.0

        return genotype, phred_qual, phred_gq

    @staticmethod
    def get_proper_alleles(positional_record, genotype):
        alts = [(positional_record.alt1, positional_record.alt1_type),
                (positional_record.alt2, positional_record.alt2_type)]

        gts = genotype.split('/')
        refined_alt = []

        if gts[0] == '0' and gts[1] == '0':
            refined_alt.append('.')
        if gts[0] == '1' or gts[1] == '1':
            refined_alt.append(alts[0])
        if gts[0] == '2' or gts[1] == '2':
            if len(alts) > 1:
                refined_alt.append(alts[1])
            elif genotype == '0/2':
                refined_alt.append(alts[0])
                genotype = '0/1'
            elif genotype == '2/2':
                refined_alt.append(alts[0])
                genotype = '1/1'
            elif genotype == '1/2':
                genotype = '0/1'

        if len(refined_alt) == 1:
            ref, alt1, alt2 = VCFWriter.solve_single_alt(refined_alt, positional_record.ref)
        else:
            ref, alt1, alt2 = VCFWriter.solve_multiple_alts(refined_alt, positional_record.ref)

        refined_alt = [alt1, alt2]
        refined_gt = genotype
        if genotype == '0/2':
            refined_gt = '0/1'
        if genotype == '2/2':
            refined_gt = '1/1'

        return ref, refined_alt, refined_gt

    @staticmethod
    def get_filter(record, last_end):
        chrm, st_pos, end_pos, ref, alt_field, genotype, phred_qual, phred_gq = record
        if st_pos < last_end:
            return 'conflictPos'
        if genotype == '0/0':
            return 'refCall'
        if phred_qual < 0:
            return 'lowQUAL'
        if phred_gq < 0:
            return 'lowGQ'
        return 'PASS'

    def get_vcf_header(self, sample_name):
        header = VariantHeader()
        items = [('ID', "PASS"),
                 ('Description', "All filters passed")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "refCall"),
                 ('Description', "Call is homozygous")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowGQ"),
                 ('Description', "Low genotype quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "lowQUAL"),
                 ('Description', "Low variant call quality")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "conflictPos"),
                 ('Description', "Overlapping record")]
        header.add_meta(key='FILTER', items=items)
        items = [('ID', "GT"),
                 ('Number', 1),
                 ('Type', 'String'),
                 ('Description', "Genotype")]
        header.add_meta(key='FORMAT', items=items)
        items = [('ID', "GQ"),
                 ('Number', 1),
                 ('Type', 'Float'),
                 ('Description', "Genotype Quality")]
        header.add_meta(key='FORMAT', items=items)
        bam_sqs = self.bam_handler.get_header_sq()
        for sq in bam_sqs:
            id = sq['SN']
            ln = sq['LN']
            items = [('ID', id),
                     ('length', ln)]
            header.add_meta(key='contig', items=items)

        header.add_sample(sample_name)

        return header
