import sys
from scipy import misc

"""
This script creates pileup images given a vcf record, bam alignment file and reference fasta file.

imageChannels: Handles how many channels to create for each base and their structure

"""
MAP_QUALITY_FILTER = 5
MAX_COLOR_VALUE = 254
BASE_QUALITY_CAP = 40
MAP_QUALITY_CAP = 60
MIN_DELETE_QUALITY = 20
MATCH_CIGAR_CODE = 0
INSERT_CIGAR_CODE = 1
DELETE_CIGAR_CODE = 2

HOM_CLASS = 0
HET_CLASS = 1
HOM_ALT_CLASS = 2


class ImageGenerator:
    """
    Processes a pileup around a position
    """
    @staticmethod
    def get_class_label_for_alt1(gt):
        h1, h2 = gt
        if h1 == 1 or h2 == 1:
            if h1 == h2:
                return HOM_ALT_CLASS
            else:
                return HET_CLASS
        return HOM_CLASS

    @staticmethod
    def get_class_label_for_alt2(gt):
        h1, h2 = gt
        if h1 == 2 or h2 == 2:
            if h1 == h2:
                return HOM_ALT_CLASS
            else:
                return HET_CLASS
        return HOM_CLASS

    @staticmethod
    def get_class_label_for_combined_alt(gt):
        h1, h2 = gt
        if h1 == 0 and h2 == 0:
            return HOM_CLASS

        if h1 == 0 or h2 == 0:
            return HET_CLASS

        return HOM_ALT_CLASS

    @staticmethod
    def get_combined_records_for_two_alts(record):
        """
        Returns records for sites where we have two alternate alleles.
        :param record: Record that belong to the site
        :return: Records of a site
        """
        chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type_alt1, rec_type_alt2 = record[0:8]
        # get the genotypes from the record
        gt = record[-1]
        gt1 = ImageGenerator.get_class_label_for_alt1(gt)
        gt2 = ImageGenerator.get_class_label_for_alt2(gt)
        # get the genotype of the images where both of these alleles are used together
        gt3 = ImageGenerator.get_class_label_for_combined_alt(gt)

        # create two separate records for each of the alleles
        rec_1 = [chr_name, pos_start, pos_end, ref, alt1, '.', rec_type_alt1, 0, gt1]
        rec_2 = [chr_name, pos_start, pos_end, ref, alt2, '.', rec_type_alt2, 0, gt2]
        if gt3 is not None:
            # if gt3 is not invalid create the record where both of the alleles are used together
            rec_3 = [chr_name, pos_start, pos_end, ref, alt1, alt2, rec_type_alt1, rec_type_alt2, gt3]
            return [rec_1, rec_2, rec_3]

        return [rec_1, rec_2]

    @staticmethod
    def save_image_as_png(pileup_array, save_dir, file_name):
        pileup_array_2d = pileup_array.reshape((pileup_array.shape[0], -1))
        try:
            misc.imsave(save_dir + file_name + ".png", pileup_array_2d, format="PNG")
        except:
            sys.stderr.write("ERROR SAVING FILE: " + save_dir + file_name + ".png" + "\n")

    @staticmethod
    def generate_and_save_candidate_images(candidate_list, summary_writer):
        # declare the size of the image
        if len(candidate_list) == 0:
            return

        # list of image records to be generated
        image_record_set = []
        # expand the records for sites where two alleles are found
        for candidate in candidate_list:
            image_record_set.extend(candidate.get_candidate_record())

        # index of the image we generate the images
        indx = 0
        for img_record in image_record_set:
            summary_writer.write(img_record + '\n')
            indx += 1
