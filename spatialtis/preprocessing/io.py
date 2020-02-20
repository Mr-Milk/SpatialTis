"""Convert images and masks to anndata object
"""

import re
from pathlib import Path
from typing import Optional, Sequence, Union

import anndata as ad
import numpy as np
import pandas as pd
import ray
from skimage.external import tifffile
from skimage.io import imread

from spatialtis.config import ISOTOPES_MASS_NUMBER_MAP, ISOTOPES_NAME

from ._utils import (config, config_file, filter_channels, get_cell_exp_stack,
                     mask2cells)
from .geom import geom_cells


class read_all_ROIs:
    """read all rois, one roi one folder

    Organize your folder as your experiment conditions.
    We will exhaust the directory.
    Length of condition_name must be the same as the depth of your folders.

    Args:
        entry: the entry of your folders
        conditions: names for each level of your folders
        mask_pattern: name pattern for mask image

    Example:
        Let's say your file structure look like this::

            Data
            ├── Patient1
            │   ├── Sample1
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   ├──Sample2
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   └──Sample3
            │       ├── ROI1
            │       └── ROI2
            ├── Patient2
            ├── Patient3
            └── ...


    Read all your data::

        entry = '.\\Data'
        condition_name = ['Patient', 'Sample', 'ROI']
        metadata = '.\\metadata.csv'

        data = read_all_ROIs(entry, condition_name)
                .config_file(metadata, channel_col='channels', marker_col='markers')

        # set mp=True for parallelism
        data = data.to_anndata(mp=True)

    """

    def __init__(
        self,
        entry: Union[Path, str],
        conditions: Sequence,
        mask_pattern: str = "*mask*",
    ):
        self.__mask_pattern = mask_pattern
        self.__channels = list()
        self.__markers = dict()
        self.__channels_files = dict()

        self.__tree = []
        self.__obs_name = conditions
        self.__depth = len(conditions)
        self.__exhaust_dir(entry)
        self.__obs = [p.parts[-self.__depth:] for p in self.__tree]
        self.__var = None

    # walk through the directory, until there is no directory
    def __exhaust_dir(
        self, path: Union[Path, str],
    ):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self.__tree.append(f)
            if f.parent in self.__tree:
                self.__tree.remove(f.parent)
            self.__exhaust_dir(f)

    @property
    def tree(self):
        return self.__obs

    @property
    def vars(self):
        return self.__var

    def config_file(
        self,
        metadata: Union[Path, str],
        channel_col: Optional[str] = None,
        marker_col: Optional[str] = None,
        sep: str = ",",
    ):
        """config with file

        A config file should contain at least one column that tells channels.

        Args:
            metadata: path to your config file
            channel_col: column name of channels
            marker_col: column name of markers
            sep: the delimiter of your file, eg.',' for `.csv`, ' ' for `.txt`, ...

        Returns:
            self

        """
        # we use callback function to set modify some info after config
        return config_file(
            self,
            metadata,
            channel_col=channel_col,
            marker_col=marker_col,
            sep=sep,
            callback=set_info,
        )

    def config(self, channels=None, markers=None):
        # we use callback function to set modify some info after config
        return config(self, channels=channels, markers=markers, callback=set_info)

    @property
    def channels(self):
        return self.__channels

    @property
    def markers(self):
        return self.__markers

    def to_anndata(self, polygonize=True, method="mean", mp=False):

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
            for i, d in enumerate(self.__tree):
                results.append(
                    _get_roi.remote(
                        d,
                        self.__channels,
                        self.__markers,
                        polygonize,
                        method,
                        self.__obs[i],
                    )
                )

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
            for i, d in enumerate(self.__tree):

                if len(self.__markers) >= 1:
                    roi = read_ROI(d).config(
                        channels=self.__channels, markers=self.__markers
                    )
                else:
                    roi = read_ROI(d).config(channels=self.__channels)

                exp, cells = roi.exp_matrix(polygonize=polygonize, method=method)

                cell_count = len(cells[0])
                obs = np.repeat(np.array([self.__obs[i]]), cell_count, axis=0)

                X += exp
                ann_obs += list(obs)
                areas += cells[0]
                shapes += cells[1]
                centroids += cells[2]
                eccentricities += cells[3]

        # print(len(ann_obs), len(areas))
        # anndata require str index, hard set to str
        ann_obs = pd.DataFrame(
            ann_obs,
            columns=self.__obs_name,
            index=[str(i) for i in range(0, len(ann_obs))],
        )
        ann_obs["area"] = areas
        ann_obs["cell_shape"] = shapes
        ann_obs["centroid"] = centroids
        ann_obs["eccentricity"] = eccentricities

        X = np.asarray(X, dtype=float)

        return ad.AnnData(X, obs=ann_obs, var=self.__var, dtype="float")


def set_info(cls):
    lc = len(cls._channels)
    lm = len(cls._markers)

    if lc == 0:
        try:
            cls._channels = read_ROI(cls._tree[0]).channels
        finally:
            cls._var = pd.DataFrame({"Channels": cls._channels})
    elif (lc > 0) & (lm > 0):
        cls._var = pd.DataFrame(
            {"Channels": cls._channels, "Markers": list(cls._markers.values())}
        )
    elif (lc > 0) & (lm == 0):
        cls._var = pd.DataFrame({"Channels": cls._channels})
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
        self.__mask_img = mask[0]

        # get channels info
        self.__channels = list()
        self.__markers = dict()
        self.__channels_files = dict()
        for img in Path(folder).iterdir():
            if img not in mask:
                with tifffile.TiffFile(str(img)) as tif:
                    col_name = tif.pages[0].tags["page_name"].value.decode()
                    pattern = re.compile(r"([a-zA-Z]+)([0-9]{2,})")
                    col = re.findall(pattern, col_name)
                    correct_isotopes = False
                    for c in col:
                        if c[0] in ISOTOPES_NAME:
                            if int(c[1]) in ISOTOPES_MASS_NUMBER_MAP[c[0]]:
                                cname = c[0] + c[1]
                                self.__channels.append(cname)
                                self.__channels_files[cname] = img
                                correct_isotopes = True
                                break
                    if not correct_isotopes:
                        print(f"Your channel isotope not exists. File: {img}")

    def config_file(self, metadata, channel_col=None, marker_col=None, sep=","):
        return config_file(
            self, metadata, channel_col=channel_col, marker_col=marker_col, sep=sep
        )

    def config(self, channels=None, markers=None):
        selected_channels = filter_channels(self, channels=channels)
        return config(self, channels=selected_channels, markers=markers)

    @property
    def channels(self):
        return self.__channels

    @property
    def markers(self):
        return self.__markers

    def exp_matrix(self, polygonize=True, method="mean"):
        if len(self.__markers) == 0:
            print("ATTENTION: NO marker specific, using channels' name instead.")
        cells = mask2cells(self.__mask_img)

        stacks = np.asarray([imread(self.__channels_files[c]) for c in self.__channels])
        data = get_cell_exp_stack(stacks, cells, method=method)
        print(f"Detected {len(data)} cells.")
        if polygonize:
            return data, geom_cells(cells)
        else:
            return data, cells
