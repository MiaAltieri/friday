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
        # "Some images referenced in the CSV file were not found"
        self.transform = transform

        self.file_info = list(data_frame[0])
        self.dict_info = list(data_frame[1])
        self.record = list(data_frame[2])

    def __getitem__(self, index):
        # load the image
        label_file_path = self.file_info[index] + ".label"
        # load the labels
        label = torch.load(label_file_path)
        label = label.type(torch.LongTensor)

        dict_path = self.dict_info[index]
        record = self.record[index]

        return label, record, dict_path

    def __len__(self):
        return len(self.file_info)
