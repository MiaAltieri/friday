import h5py
import yaml
import numpy as np
from modules.python.Options import ImageSizeOptions


class DataStore(object):
    """Class to read/write to a FRIDAY's file"""
    _prediction_path_ = 'predictions'
    _groups_ = 'friday_predictions'

    def __init__(self, filename, mode='r'):
        self.filename = filename
        self.mode = mode

        self._sample_keys = set()
        self.file_handler = h5py.File(self.filename, self.mode)

        self._meta = None

    def __exit__(self, *args):
        if self.mode != 'r' and self._meta is not None:
            self._write_metadata(self.meta)
        self.file_handler.close()

    def _write_metadata(self, data):
        """Save a data structure to file within a yml str."""
        for group, d in data.items():
            if group in self.file_handler:
                del self.file_handler[group]
            self.file_handler[group] = yaml.dump(d)

    def _load_metadata(self, groups=None):
        """Load meta data"""
        if groups is None:
            groups = self._groups_
        return {g: yaml.load(self.file_handler[g][()]) for g in groups if g in self.file_handler}

    @property
    def meta(self):
        if self._meta is None:
            self._meta = self._load_metadata()
        return self._meta

    def update_meta(self, meta):
        """Update metadata"""
        self._meta = self.meta
        self._meta.update(meta)

    def write_predictions(self, chromosome_names, image_names, predictions):
        if 'friday_predictions' not in self.meta:
            self.meta['friday_predictions'] = set()

        for image_name, chromosome_name, prediction in zip(image_names, chromosome_names, predictions):
            if image_name not in self.meta['friday_predictions']:
                # update the meta
                self.meta['friday_predictions'].add(image_name)
                # save all attributes
                self.file_handler['{}/{}/{}'.format(self._prediction_path_, chromosome_name, image_name)] = prediction
