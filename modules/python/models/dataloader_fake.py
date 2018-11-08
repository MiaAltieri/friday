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
        self.file_index = list(data_frame[1])
        self.dict_info = list(data_frame[2])
        self.record = list(data_frame[3])

    def __getitem__(self, index):
        dict_path = self.dict_info[index]
        record = self.record[index]

        hdf5_image = self.file_info[index]
        hdf5_index = int(self.file_index[index])
        hdf5_file = h5py.File(hdf5_image, 'r')

        label_dataset = hdf5_file['labels']
        label = np.array(label_dataset[hdf5_index], dtype=np.long)
        label = torch.from_numpy(label).type(torch.LongTensor)

        return label, record, dict_path

    def __len__(self):
        return len(self.file_info)
