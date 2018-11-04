import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import h5py
import torch
import pickle


class SequenceDataset(Dataset):
    """
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        data_frame = pd.read_csv(csv_path, header=None, dtype=str)
        # assert data_frame[0].apply(lambda x: os.path.isfile(x.split(' ')[0])).all(), \
        #     "Some images referenced in the CSV file were not found"
        self.transform = transform
        self.file_info = list(data_frame[0])
        self.dict_info = list(data_frame[1])
        self.record = list(data_frame[2])

    @staticmethod
    def load_dictionary(dictionary_location):
        f = open(dictionary_location, 'rb')
        dict = pickle.load(f)
        f.close()
        return dict

    def __getitem__(self, index):
        # load the image
        image_file_path = self.file_info[index] + ".image"
        label_file_path = self.file_info[index] + ".label"

        image = torch.load(image_file_path)
        image = image.transpose(1, 2)
        # load the labels
        label = torch.load(label_file_path)

        image = image.type(torch.FloatTensor)
        label = label.type(torch.LongTensor)

        return image, label

    def __len__(self):
        return len(self.file_info)
