import h5py
import yaml
import numpy as np
from modules.python.Options import ImageSizeOptions


class DataStore(object):
    """Class to read/write to a FRIDAY's file"""
    _image_path_ = 'images'
    _candidate_path_ = 'candidates'
    _prediction_path_ = 'predictions'
    _groups_ = ('friday_candidates', 'friday_images', 'friday_predictions')

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

    def write_images(self, images, chromosome_name):
        if 'friday_images' not in self.meta:
            self.meta['friday_images'] = set()

        image_datatype = np.dtype([('chromosome_name', h5py.special_dtype(vlen=str)),
                                   ('start_pos', np.uint64),
                                   ('end_pos', np.uint64),
                                   ('label', np.uint8),
                                   ('name', h5py.special_dtype(vlen=str)),
                                   ('image', h5py.special_dtype(vlen=np.uint8))])

        data_array = np.array(images, dtype=image_datatype)
        dataset = self.file_handler.create_dataset('{}/{}/'.format(self._image_path_, chromosome_name),
                                                   (len(images),), dtype=image_datatype)

        dataset[...] = data_array

    def write_candidates(self, candidates, chromosome_name):
        if 'friday_candidates' not in self.meta:
            self.meta['friday_candidates'] = set()

        candidate_datatype = np.dtype([('chromosome_name', h5py.special_dtype(vlen=str)),
                                       ('pos_start', np.uint64),
                                       ('pos_end', np.uint64),
                                       ('name', h5py.special_dtype(vlen=str)),
                                       ('ref', h5py.special_dtype(vlen=str)),
                                       ('alternate_alleles', h5py.special_dtype(vlen=str)),
                                       ('allele_depths', h5py.special_dtype(vlen=np.uint32)),
                                       ('allele_frequencies', h5py.special_dtype(vlen=np.float32)),
                                       ('genotype', h5py.special_dtype(vlen=np.uint8)),
                                       ('image_names', h5py.special_dtype(vlen=str))])

        data_array = np.array(candidates, dtype=candidate_datatype)
        dataset = self.file_handler.create_dataset('{}/{}/'.format(self._candidate_path_, chromosome_name),
                                                   (len(candidates),), dtype=candidate_datatype)

        dataset[...] = data_array

