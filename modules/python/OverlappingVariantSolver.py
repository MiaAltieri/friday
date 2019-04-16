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

    def solve_variant_group(self, grouped_variant):
        if len(grouped_variant) == 1:
            return grouped_variant[0]
        else:
            print(grouped_variant)

    def solve_overlapping_variants(self, sorted_called_variants):
        grouped_variants = self.group_variants(sorted_called_variants)
        solved_variants = []
        for group in grouped_variants:
            solved_variants.append(self.solve_variant_group(group))

