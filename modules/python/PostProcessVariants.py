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
from modules.python.OverlappingVariantSolver import OverLappingVariantSolver
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq')


class PostProcessVariants:
    @staticmethod
    def get_file_paths_from_directory(directory_path):
        """
        Returns all paths of files given a directory path
        :param directory_path: Path to the directory
        :return: A list of paths of files
        """
        file_paths = [join(directory_path, file) for file in listdir(directory_path)
                      if isfile(join(directory_path, file))
                      and file[-3:] == 'hdf']
        return file_paths

    def get_candidates(self, candidate_hdf5_directory):
        candidate_files = self.get_file_paths_from_directory(candidate_hdf5_directory)
        candidate_by_chromosome = defaultdict(list)

        for file_name in candidate_files:
            hdf5_file = h5py.File(file_name, 'r')

            chromosome_names = hdf5_file['candidates'].keys()
            for chromosome_name in chromosome_names:
                for i in range(0, hdf5_file['candidates'][chromosome_name].shape[0]):
                    candidate_by_chromosome[chromosome_name].append(hdf5_file['candidates'][chromosome_name][i])

            hdf5_file.close()

        return candidate_by_chromosome

    @staticmethod
    def get_predictions(prediction_file):
        image_name_prediction_dict = defaultdict(lambda: defaultdict(list))
        hdf5_file = h5py.File(prediction_file, 'r')

        chromosome_names = hdf5_file['predictions'].keys()
        for chromosome_name in chromosome_names:
            image_file_names = hdf5_file['predictions'][chromosome_name].keys()
            for image_name in image_file_names:
                prediction_list = np.array(hdf5_file['predictions'][chromosome_name][image_name]).tolist()
                image_name_prediction_dict[chromosome_name][image_name] = prediction_list

        hdf5_file.close()

        return image_name_prediction_dict

    @staticmethod
    def get_quals(predictions, prediction_index):
        perror = 1.0 - min(predictions[prediction_index], 1.0 - 1e-15)
        if perror < 0.0 or perror > 1.0:
            # need to raise stuff, this should be replaced
            sys.stderr.write(TextColor.RED + "ERROR: INVALID PERROR VALUE: " + str(perror) + TextColor. END)
            exit()

        gq = -10.0 * math.log10(perror)

        perrqual = 1.0 - min(sum(predictions[1:]), 1.0 - 1e-15)
        if perrqual < 0.0 or perrqual > 1.0:
            # need to raise stuff, this should be replaced
            sys.stderr.write(TextColor.RED + "ERROR: INVALID QUAL VALUE: " + str(perrqual) + TextColor. END)
            exit()
        qual = -10.0 * math.log10(perrqual)
        rounded_qual = round(qual, 8)

        return gq, rounded_qual

    def get_genotype_for_single_allelic_site(self, predictions):
        genotype_index = np.argmax(np.array(predictions))

        gt = [0, 0]
        if genotype_index == 1:
            gt = [0, 1]
        elif genotype_index == 2:
            gt = [1, 1]

        gq, qual = self.get_quals(np.array(predictions), genotype_index)
        return gt, gq, qual

    def get_callable_alleles_from_multialleleic_sites(self, predictions):
        callable_alleles = []
        removable_alleles = []
        highest_qual = -1
        highest_qual_allele_index = -1
        for allele_indices, prediction in predictions:
            if len(allele_indices) == 1:
                genotype_index = np.argmax(np.array(prediction))
                gq, qual = self.get_quals(np.array(prediction), genotype_index)

                if qual < VariantPostProcessingOptions.MULTI_ALLELE_QUAL_THRESHOLD:
                    removable_alleles.append(allele_indices[0])
                else:
                    callable_alleles.append(allele_indices[0])

                if highest_qual < 0 or qual > highest_qual:
                    highest_qual = qual
                    highest_qual_allele_index = allele_indices[0]

        if not callable_alleles:
            return [highest_qual_allele_index]

        return callable_alleles

    def get_genotype_for_multi_allelic_site(self, predictions, callable_alleles):
        allele_combination_predictions = defaultdict(lambda: [])

        for allele_indices, prediction in predictions:
            ref = [0]
            alt = [i+1 for i in allele_indices]

            valid_set = True
            for indx in allele_indices:
                if indx not in callable_alleles:
                    valid_set = False
                    break

            if not valid_set:
                continue

            p_hom, p_het, p_homalt = tuple(prediction)

            for (alt0, alt1, p) in [(ref, ref, p_hom), (ref, alt, p_het), (alt, alt, p_homalt)]:
                for alt_combination in itertools.product(alt0, alt1):
                    alt_combination = tuple(sorted(alt_combination))
                    allele_combination_predictions[alt_combination].append(p)

        genotype_predictions = defaultdict(float)
        total_prob = 0.0
        for genotype in allele_combination_predictions.keys():
            genotype_predictions[genotype] = min(allele_combination_predictions[genotype])
            total_prob += genotype_predictions[genotype]

        predictions = []
        genotypes = []
        for genotype in genotype_predictions:
            if total_prob <= 0.0:
                genotype_predictions[genotype] = 1.0 / float(len(genotype_predictions.keys()))
            else:
                genotype_predictions[genotype] = float(genotype_predictions[genotype]) / float(total_prob)

            genotypes.append(genotype)
            predictions.append(genotype_predictions[genotype])

        genotype_index = int(np.argmax(np.array(predictions)))
        gt = list(genotypes[genotype_index])
        gq, qual = self.get_quals(predictions, genotype_index)

        return gt, gq, qual

    def get_canonical_variants_from_candidates(self, candidate_set, image_name_to_prediction):

        all_called_candidates = []
        for candidate in candidate_set:
            chromosome_name, pos_start, pos_end, name, ref, alternate_alleles, allele_depths, \
                             allele_frequencies, candidate_genotype, image_names = candidate

            # make alternate allele list as a list of strings again
            alternate_alleles = alternate_alleles.strip().replace("[", "").replace("]", "").replace("'", "").split(' ')

            # create a map of allelelic prediction combination
            predictions = []
            image_names = image_names.strip().replace("[", "").replace("]", "").replace("'", "").split(' ')
            for image_name in image_names:
                if image_name.strip() not in image_name_to_prediction.keys():
                    sys.stderr.write(TextColor.RED + "ERROR: IMAGE PREDICTION DOES NOT EXIST IN DICTIONARY: " +
                                     image_name + TextColor. END)
                    exit()

                genotype_predictions = image_name_to_prediction[image_name.strip()]
                allele_indices = [int(i) for i in image_name.strip().split('_')[-1]]
                predictions.append((allele_indices, genotype_predictions))

            if len(alternate_alleles) == 1:
                # THESE ARE EASY TO SOLVE JUST CONVERT THE PREDICTION TO QUAL VALUES
                # AND SEE IF THE VALUE IS GREATER OR NOT

                genotype, gq, qual = self.get_genotype_for_single_allelic_site(predictions[0][1])
                if genotype == [0, 0] or qual < VariantPostProcessingOptions.SINGLE_QUAL_THRESHOLD:
                    continue
                else:
                    called_variant = Candidate(chromosome_name, pos_start, pos_end, ref, alternate_alleles,
                                               allele_depths.tolist(), allele_frequencies.tolist(), genotype,
                                               qual, gq)
                    all_called_candidates.append(called_variant)
            else:
                callable_alleles = self.get_callable_alleles_from_multialleleic_sites(predictions)
                genotype, gq, qual = self.get_genotype_for_multi_allelic_site(predictions, callable_alleles)
                if genotype == [0, 0] or qual < VariantPostProcessingOptions.SINGLE_QUAL_THRESHOLD:
                    continue
                else:
                    called_variant = Candidate(chromosome_name, pos_start, pos_end, ref, alternate_alleles,
                                               allele_depths.tolist(), allele_frequencies.tolist(), genotype,
                                               qual, gq)
                    all_called_candidates.append(called_variant)

        sorted_candidates = sorted(all_called_candidates, key=lambda element: (element[0], element[1], element[2]))
        return sorted_candidates

    def perform_post_processing(self, candidate_hdf5_directory, prediction_hdf_file):
        sys.stderr.write(TextColor.GREEN + "INFO: LOADING CANDIDATES FROM FILES\n" + TextColor.END)
        candidate_by_chromosome = self.get_candidates(candidate_hdf5_directory)
        sys.stderr.write(TextColor.GREEN + "INFO: LOADING PREDICTIONS FROM FILE\n" + TextColor.END)
        image_name_to_predictions = self.get_predictions(prediction_hdf_file)
        sys.stderr.write(TextColor.GREEN + "INFO: DATA LOADING COMPLETE\n" + TextColor.END)

        # this can be easily parallelized
        for chromosome_name in candidate_by_chromosome.keys():
            sys.stderr.write(TextColor.GREEN + "INFO: PROCESSING " + str(chromosome_name) + "\n" + TextColor.END)
            sorted_called_variants = self.get_canonical_variants_from_candidates(
                                                                            candidate_by_chromosome[chromosome_name],
                                                                            image_name_to_predictions[chromosome_name])

            sys.stderr.write(TextColor.GREEN + "INFO: TOTAL " + str(len(sorted_called_variants)) +
                             " VARIANTS FOUND IN " + chromosome_name + " INCLUDING OVERLAPPING\n" + TextColor.END)

            sys.stderr.write(TextColor.GREEN + "INFO: SOLVING OVERLAPPING VARIANTS\n" + TextColor.END)
            overlap_solver = OverLappingVariantSolver()
            overlap_solver.solve_overlapping_variants(sorted_called_variants)


