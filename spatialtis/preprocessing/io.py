import re
import numpy as np
import pandas as pd
import anndata as ad
from skimage.io import imread
from pathlib import Path
from skimage.external import tifffile

import ray

from typing import Sequence, Union

from .geom import geom_cells

from spatialtis.config import ISOTOPES_NAME, ISOTOPES_MASS_NUMBER_MAP
from ._utils import config, config_file, mask2cells, get_cell_exp_stack, filter_channels


class read_all_ROIs:

    def __init__(self, entry, conditions, mask_pattern="*mask*"):
        """
        Organize your folder as your experiment conditions.
        We will exhaust the directory.
        Length of condition_name must be the same as the depth of your folders.
    
        Example:

        condition_name = ['Patient', 'Sample', 'ROI']

        Entry
        |------Patient1
            |------Sample1
                |------ROI1
                |------ROI2
            |------Sample2
                |------ROI1
                |------ROI2
            |------Sample3
                |------ROI1
                |------ROI2
        |------Patient2
            |------Sample1
            |------Sample2
            |------Sample3
        |------Patient3
            |------Sample1
            |------Sample2
            |------Sample3
        
        """

        self._mask_pattern = mask_pattern
        self._channels = list()
        self._markers = dict()
        self._channels_files = dict()

        self._tree = []
        self._obs_name = conditions
        self._depth = len(conditions)
        self._exhaust_dir(entry)
        self._obs = [p.parts[-self._depth:] for p in self._tree]
        self._var = None

    # walk through the directory, until there is no directory
    def _exhaust_dir(self, path):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self._tree.append(f)
            if f.parent in self._tree:
                self._tree.remove(f.parent)
            self._exhaust_dir(f)

    @property
    def tree(self):
        return self._obs

    @property
    def vars(self):
        return self._var

    def config_file(self, metadata, channel_col=None, marker_col=None, sep=','):
        # we use callback function to set modify some info after config
        return config_file(self, metadata, channel_col=channel_col, marker_col=marker_col, sep=sep, callback=set_info)

    def config(self, channels=None, markers=None):
        # we use callback function to set modify some info after config
        return config(self, channels=channels, markers=markers, callback=set_info)

    @property
    def channels(self):
        return self._channels

    @property
    def markers(self):
        return self._markers

    def to_anndata(self, polygonize=True, method='mean', mp=False):

        X = []
        ann_obs = []

        areas = []
        shapes = []
        centroids = []
        eccentricities = []

        if mp:
            @ray.remote
            def _get_roi(t, channels, markers, pg, mt, obsi):
                if len(markers) >= 1:
                    roi = read_ROI(t).config(channels=channels, markers=markers)
                else:
                    roi = read_ROI(t).config(channels=channels)

                exp, cells = roi.exp_matrix(polygonize=pg, method=mt)

                cell_count = len(cells[0])
                obs = np.repeat(np.array([obsi]), cell_count, axis=0)
                print(f"Added: {' '.join(obsi)}")
                return [exp, list(obs), cells]

            results = []
            for i, d in enumerate(self._tree):
                results.append(_get_roi.remote(d, self._channels, self._markers, polygonize, method, self._obs[i]))

            results = ray.get(results)

            for i in results:
                X += i[0]
                ann_obs += i[1]
                cells = i[2]
                areas += cells[0]
                shapes += cells[1]
                centroids += cells[2]
                eccentricities += cells[3]

        else:
            for i, d in enumerate(self._tree):

                if len(self._markers) >= 1:
                    roi = read_ROI(d).config(channels=self._channels, markers=self._markers)
                else:
                    roi = read_ROI(d).config(channels=self._channels)

                exp, cells = roi.exp_matrix(polygonize=polygonize, method=method)

                cell_count = len(cells[0])
                obs = np.repeat(np.array([self._obs[i]]), cell_count, axis=0)

                X += exp
                ann_obs += list(obs)
                areas += cells[0]
                shapes += cells[1]
                centroids += cells[2]
                eccentricities += cells[3]

        # print(len(ann_obs), len(areas))
        # anndata require str index, hard set to str
        ann_obs = pd.DataFrame(ann_obs, columns=self._obs_name, index=[str(i) for i in range(0, len(ann_obs))])
        ann_obs['area'] = areas
        ann_obs['cell_shape'] = shapes
        ann_obs['centroid'] = centroids
        ann_obs['eccentricity'] = eccentricities

        X = np.asarray(X, dtype=float)

        return ad.AnnData(X, obs=ann_obs, var=self._var, dtype='float')


def set_info(cls):
    lc = len(cls._channels)
    lm = len(cls._markers)

    if lc == 0:
        try:
            cls._channels = read_ROI(cls._tree[0]).channels
        finally:
            cls._var = pd.DataFrame({'Channels': cls._channels})
    elif (lc > 0) & (lm > 0):
        cls._var = pd.DataFrame({'Channels': cls._channels, 'Markers': list(cls._markers.values())})
    elif (lc > 0) & (lm == 0):
        cls._var = pd.DataFrame({'Channels': cls._channels})
    # anndata require str index, hard set everything to str
    cls._var.index = [str(i) for i in range(0, len(cls._channels))]


class read_ROI:
    """
    Your .tif/.tiff file should be exported from MCD viewer
    or you can specific channel name in 'page_name' field.
    """

    def __init__(self, folder, mask_pattern="*mask*"):
        """
        Specific the mask image, or automatically select the img name contain "mask".
        """
        self._work_dir = Path(folder)

        # try to find mask
        mask = [i for i in self._work_dir.glob(mask_pattern)]
        if len(mask) == 0:
            raise FileNotFoundError("No mask found")
        elif len(mask) > 1:
            print(f"Found more than one mask, use {mask[0].name}")
        self._mask_img = mask[0]

        # get channels info
        self._channels = list()
        self._markers = dict()
        self._channels_files = dict()
        for img in Path(folder).iterdir():
            if img not in mask:
                with tifffile.TiffFile(str(img)) as tif:
                    col_name = tif.pages[0].tags['page_name'].value.decode()
                    pattern = re.compile(r'([a-zA-Z]+)([0-9]{2,})')
                    col = re.findall(pattern, col_name)
                    correct_isotopes = False
                    for c in col:
                        if c[0] in ISOTOPES_NAME:
                            if int(c[1]) in ISOTOPES_MASS_NUMBER_MAP[c[0]]:
                                cname = c[0] + c[1]
                                self._channels.append(cname)
                                self._channels_files[cname] = img
                                correct_isotopes = True
                                break
                    if not correct_isotopes:
                        print(f"Your channel isotope not exists. File: {img}")

    def config_file(self, metadata, channel_col=None, marker_col=None, sep=','):
        return config_file(self, metadata, channel_col=channel_col, marker_col=marker_col, sep=sep)

    def config(self, channels=None, markers=None):
        selected_channels = filter_channels(self, channels=channels)
        return config(self, channels=selected_channels, markers=markers)

    @property
    def channels(self):
        return self._channels

    @property
    def markers(self):
        return self._markers

    def exp_matrix(self, polygonize=True, method='mean'):
        if len(self._markers) == 0:
            print("ATTENTION: NO marker specific, using channels' name instead.")
        cells = mask2cells(self._mask_img)

        stacks = np.asarray([imread(self._channels_files[c]) for c in self._channels])
        data = get_cell_exp_stack(stacks, cells, method=method)
        print(f"Detected {len(data)} cells.")
        if polygonize:
            return data, geom_cells(cells)
        else:
            return data, cells
