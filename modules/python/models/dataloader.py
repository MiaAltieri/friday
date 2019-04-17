import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torchvision.transforms as transforms
from os.path import isfile, join
from os import listdir
import h5py
import torch
import pickle


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-3:] == 'hdf']
    return file_paths


class SequenceDataset(Dataset):
    """
    Arguments: A directory to all the images
    """

    def __init__(self, image_directory, transform=None):
        hdf_files = get_file_paths_from_directory(image_directory)
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []
        for hdf5_filepath in hdf_files:
            hdf5_file = h5py.File(hdf5_filepath, 'r')
            chromosome_names = hdf5_file['images'].keys()
            for chromosome_name in chromosome_names:
                for i in range(0, hdf5_file['images'][chromosome_name]['images'].shape[0]):
                    file_image_pair.append((hdf5_filepath, chromosome_name, i))
            hdf5_file.close()

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, chromosome_name, hdf5_index = self.all_images[index]

        hdf5_file = h5py.File(hdf5_filepath, 'r')
        image = hdf5_file['images'][chromosome_name]['images'][hdf5_index]

        image = np.array(image, dtype=np.uint8)
        label = hdf5_file['images'][chromosome_name]['labels'][hdf5_index]
        label = np.array(label, dtype=np.long)
        hdf5_file.close()

        image = self.transform(image)
        image = image.transpose(1, 2)

        label = torch.from_numpy(label).type(torch.LongTensor)

        return image, label

    def __len__(self):
        return len(self.all_images)
