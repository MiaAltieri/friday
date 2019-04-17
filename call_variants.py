import argparse
import sys
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader
from modules.python.models.ModelHander import ModelHandler
from modules.python.models.dataloader_predict import SequenceDataset
from modules.python.Options import ImageSizeOptions
from modules.python.TextColor import TextColor
from modules.python.DataStore_Prediction import DataStore
from tqdm import tqdm
import os
"""
This script uses a trained model to call variants on a given set of images generated from the genome.
The process is:
- Create a prediction table/dictionary using a trained neural network
- Convert those predictions to a VCF file

INPUT:
- A trained model
- Set of images for prediction

Output:
- A HDF5 file containing all the variant predictions with image name.
"""


def call_variants(hdf5_filepath, batch_size, model_path, gpu_mode, num_workers, output_path):
    """
    Create a prediction table/dictionary of an images set using a trained model.
    :param hdf5_filepath: File to predict on
    :param batch_size: Batch size used for prediction
    :param model_path: Path to a trained model
    :param gpu_mode: If true, predictions will be done over GPU
    :param num_workers: Number of workers to be used by the dataloader
    :return: Prediction dictionary
    """
    # the prediction table/dictionary
    sys.stderr.write(TextColor.PURPLE + 'Loading data\n' + TextColor.END)

    test_dset = SequenceDataset(hdf5_filepath)
    testloader = DataLoader(test_dset,
                            batch_size=batch_size,
                            shuffle=False,
                            num_workers=num_workers
                            )

    sys.stderr.write(TextColor.PURPLE + 'Data loading finished\n' + TextColor.END)

    # load the model
    model, prev_ite = ModelHandler.load_model_for_training(model_path,
                                                           input_channels=ImageSizeOptions.IMAGE_CHANNELS,
                                                           num_classes=ImageSizeOptions.TOTAL_LABELS)

    if gpu_mode:
        model = torch.nn.DataParallel(model).cuda()

    # Change model to 'eval' mode (BN uses moving mean/var).
    model.eval()

    data_file_name = output_path + "friday" + "_predictions" + ".hdf"
    data_file = DataStore(data_file_name, mode='w')

    sys.stderr.write(TextColor.PURPLE + 'MODEL LOADED\n' + TextColor.END)
    with torch.no_grad():
        for images, image_names, chromosome_names in tqdm(testloader, ncols=50):
            if gpu_mode:
                images = images.cuda()

            output_ = model(images)
            m = nn.Softmax(dim=1)
            soft_probs = m(output_)
            output_preds = soft_probs.cpu()
            predictions = np.array(output_preds).tolist()
            data_file.write_predictions(chromosome_names, image_names, predictions)


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
        help="CSV directory containing all image segments for prediction."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        required=False,
        default=100,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        required=False,
        default=4,
        help="Batch size for testing, default is 100."
    )
    parser.add_argument(
        "--model_path",
        type=str,
        default='./CNN.pkl',
        help="Saved model path."
    )
    parser.add_argument(
        "--gpu_mode",
        type=bool,
        default=False,
        help="If true then cuda is on."
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

    call_variants(FLAGS.image_dir,
                  FLAGS.batch_size,
                  FLAGS.model_path,
                  FLAGS.gpu_mode,
                  FLAGS.num_workers,
                  FLAGS.output_dir)
