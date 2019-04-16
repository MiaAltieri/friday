import h5py
import numpy as np
import sys
from modules.python.TextColor import TextColor
from os.path import isfile, join
from os import listdir
from collections import defaultdict
import collections
import math
import itertools
from modules.python.Options import VariantPostProcessingOptions

"""
https://github.com/google/deepvariant/blob/master/deepvariant/haplotypes.py
IMPLMENTATION EDITED FROM DEEPVARIANT'S POST-PROCESSING PIPELINE
"""


class VariantCompatibilityCalculator(object):
    def __init__(self, overlapping_variants):
        min_start = min(v.pos_start for v in overlapping_variants)
        self.variant_indices = [
            (v.pos_start - min_start, v.pos_end - min_start) for v in overlapping_variants
        ]
        self.size = max(v.pos_end - min_start for v in overlapping_variants)

    def all_variants_compatible(self, nonref_genotype_counts, ploidy=2):
        if len(nonref_genotype_counts) != len(self.variant_indices):
            raise ValueError(
                'Variant counts must have same length as variant indices.')
        if not all(0 <= cnt <= ploidy for cnt in nonref_genotype_counts):
            raise ValueError('Invalid variant allele count for ploidy {}: {}'.format(
                ploidy, nonref_genotype_counts))

        alts_in_span = np.zeros(self.size, dtype=int)
        for cnt, (start, end) in zip(nonref_genotype_counts, self.variant_indices):
            alts_in_span[start:end] += cnt
        return np.all(alts_in_span <= ploidy)


class OverLappingVariantSolver:
    @staticmethod
    def group_variants(sorted_called_variants):
        all_groups = []
        running_group = []
        prev_chromosome = None
        prev_max_end = -1
        for variant in sorted_called_variants:
            if variant.chromosome_name != prev_chromosome or variant.pos_start >= prev_max_end:
                if running_group:
                    all_groups.append(running_group)

                running_group = [variant]
                prev_chromosome = variant.chromosome_name
                prev_max_end = variant.pos_end
            else:
                running_group.append(variant)
                prev_max_end = max(prev_max_end, variant.pos_end)

        all_groups.append(running_group)

        return all_groups

    def genotype_count(self, variant):
        return sum(g > 0 for g in variant.genotype)

    def solve_variant_group(self, grouped_variant):
        if len(grouped_variant) == 1:
            return grouped_variant[0]

        calculator = VariantCompatibilityCalculator(grouped_variant)
        nonref_counts = [self.genotype_count(v) for v in grouped_variant]

        if calculator.all_variants_compatible(nonref_counts):
            return grouped_variant

        print("NOT COMPATIBLE CALLS: ")
        for variant in grouped_variant:
            print(variant)
        print("-----------------------------")


    def solve_overlapping_variants(self, sorted_called_variants):
        grouped_variants = self.group_variants(sorted_called_variants)
        solved_variants = []
        for group in grouped_variants:
            solved_variants.append(self.solve_variant_group(group))

