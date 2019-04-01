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

    def write_images(self, images):
        if 'friday_images' not in self.meta:
            self.meta['friday_images'] = set()

        for image in images:
            if image.name not in self.meta['friday_images']:
                # update the meta
                self.meta['friday_images'].add(image.name)
                # save all attributes
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'chromosome_name')] = image.chromosome_name
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'start_pos')] = image.start_pos
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'end_pos')] = image.end_pos
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'label')] = image.label
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'name')] = image.name
                self.file_handler['{}/{}/{}'.format(self._image_path_, image.name, 'image')] = image.image
                # img_dset = self.file_handler.create_dataset('{}/{}/{}'.format(self._image_path_, image.name, 'image'),
                #                                             (ImageSizeOptions.IMAGE_HEIGHT,
                #                                              ImageSizeOptions.IMAGE_LENGTH,
                #                                              ImageSizeOptions.IMAGE_CHANNELS),
                #                                              np.uint8, compression='gzip')
                # img_dset[...] = image.image

    def write_candidates(self, candidates):
        if 'friday_candidates' not in self.meta:
            self.meta['friday_candidates'] = set()
        for candidate in candidates:
            if candidate.name not in self.meta['friday_candidates']:
                # update the meta
                self.meta['friday_candidates'].add(candidate.name)
                # save all attributes
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'chromosome_name')] = np.string_(candidate.chromosome_name)
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'pos_start')] = candidate.pos_start
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'pos_end')] = candidate.pos_end
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'name')] = np.string_(candidate.name)
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'ref')] = np.string_(candidate.ref)
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'alternate_alleles')] = [np.string_(i) for i in candidate.alternate_alleles]
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'allele_depths')] = candidate.allele_depths
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'allele_frequencies')] = candidate.allele_frequencies
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'genotype')] = candidate.genotype
                self.file_handler['{}/{}/{}'.format(self._candidate_path_, candidate.name, 'image_names')] = [np.string_(i) for i in candidate.image_names]

