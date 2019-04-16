import argparse
from modules.python.PostProcessVariants import PostProcessVariants
import os


def post_process_variants(hdf5_filepath, prediction_file, output_dir):
    # the prediction table/dictionary
    post_processor = PostProcessVariants()
    post_processor.perform_post_processing(hdf5_filepath, prediction_file)


def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--image_dir",
        type=str,
        required=True,
        help="Directory containing generated image files."
    )
    parser.add_argument(
        "--prediction_file",
        type=str,
        required=True,
        help="Prediction file."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        default='vcf_output',
        help="Output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()
    FLAGS.output_dir = handle_output_directory(FLAGS.output_dir)

    post_process_variants(FLAGS.image_dir,
                          FLAGS.prediction_file,
                          FLAGS.output_dir)
