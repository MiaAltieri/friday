import os
import numpy as np
import h5py
import pandas as pd
from torch.utils.data import Dataset
import sys


class TextColor:
    """
    Defines color codes for text used to give different mode of errors.
    """
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class PileupDataset(Dataset):
    """
    Creates a pile up image from hdf5 file as specified by csv file.
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        tmp_df = pd.read_csv(csv_path, header=None)
        # assert tmp_df[0].apply(lambda x: os.path.isfile(x)).all(), \
        #     "Some images referenced in the CSV file were not found: "
        self.rec = tmp_df[0]

    def __getitem__(self, index):
        rec = self.rec[index]
        label = int(rec[-1])
        return label, rec

    def __len__(self):
        return len(self.rec.index)