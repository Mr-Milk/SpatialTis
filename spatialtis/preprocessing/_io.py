"""Convert images and masks to anndata object
"""

import re
from pathlib import Path
from typing import Optional, Sequence, Union
import warnings

import anndata as ad
import numpy as np
import pandas as pd

from skimage.external import tifffile
from skimage.io import imread

from spatialtis.config import ISOTOPES_MASS_NUMBER_MAP, ISOTOPES_NAME

from ._utils import (config, config_file, filter_channels, get_cell_exp_stack,
                     mask2cells, set_info)
from ._geom import geom_cells


class read_ROIs:
    """read all rois, one roi one folder

    Organize your folder as your experiment conditions.
    We will exhaust the directory.
    Length of condition_name must be the same as the depth of your folders.

    Args:
        entry: the entry of your folders
        conditions: names for each level of your folders
        mask_pattern: name pattern for mask image

    Attributes:
        channels:
        markers:
        obs:
        var:

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
            stacked: bool = False,
    ):
        self.channels = list()
        self.markers = dict()

        self.__stacked = stacked,
        self.__mask_pattern = mask_pattern
        self.__channels_files = dict()

        self.tree = []
        self.__obs_name = conditions
        self.__depth = len(conditions)

        self.__exhaust_dir(entry)
        self.obs = [p.parts[-self.__depth:] for p in self.tree]

        self.var = None

    # walk through the directory, until there is no directory
    def __exhaust_dir(
            self,
            path: Union[Path, str],
    ):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self.tree.append(f)
            if f.parent in self.tree:
                self.tree.remove(f.parent)
            self.__exhaust_dir(f)

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

    def to_anndata(self, method="mean", polygonize="convex", alpha=0, mp=False):

        X = []
        ann_obs = []

        areas = []
        shapes = []
        centroids = []
        eccentricities = []

        if polygonize == "concave":
            if alpha == 0:
                raise ValueError("If using concave method, alpha value should > 0, please set argument `alpha`")
            warnings.warn("Running alphashape is very slow", RuntimeWarning)
        elif polygonize == "convex":
            pass
        else:
            raise ValueError("Polygonize options are 'convex' or 'concave'")

        if mp:
            try:
                import ray
            except ImportError:
                raise ImportError("You don't have ray installed or your OS don't support ray. Please use `mp=False`")

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
            for i, d in enumerate(self.tree):
                results.append(
                    _get_roi.remote(
                        d,
                        self.channels,
                        self.markers,
                        polygonize,
                        method,
                        self.obs[i],
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
            for i, d in enumerate(self.tree):

                if len(self.markers) >= 1:
                    roi = read_ROI(d).config(
                        channels=self.channels, markers=self.markers
                    )
                else:
                    roi = read_ROI(d).config(channels=self.channels)

                exp, cells = roi.exp_matrix(polygonize=polygonize, method=method)

                cell_count = len(cells[0])
                obs = np.repeat(np.array([self.obs[i]]), cell_count, axis=0)

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

        return ad.AnnData(X, obs=ann_obs, var=self.var, dtype="float")


class read_ROI:
    """
    Your .tif/.tiff file should be exported from MCD viewer
    or you can specific channel name in 'page_name' field.
    """

    def __init__(self, folder, mask_pattern="*mask*", stacked=False):
        """
        Specific the mask image, or automatically select the img name contain "mask".
        """
        self._work_dir = Path(folder)

        # try to find mask
        mask = [i for i in self._work_dir.glob(mask_pattern)]
        if len(mask) == 0:
            raise FileNotFoundError(f"Mask not found in {str(folder)}")
        elif len(mask) > 1:
            raise ValueError(f"Found more than one mask in {str(folder)}")
        self.__mask_img = mask[0]

        # get channels info
        self.channels = list()
        self.markers = dict()
        self.__channels_files = dict()
        self.__stacks = 0

        if stacked:
            stacks_count = 0
            for img in Path(folder).iterdir():
                if img not in mask:
                    self.__stacks = imread(str(img))
                    stacks_count += 1

            if stacks_count == 0:
                raise FileNotFoundError(f"ROI image not found in {str(folder)}")
            elif stacks_count > 1:
                raise ValueError(f"Found more than one ROI image in {str(folder)}")

        else:
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
                                    self.channels.append(cname)
                                    self.__channels_files[cname] = img
                                    correct_isotopes = True
                                    break
                        if not correct_isotopes:
                            print(f"Your channel isotope not exists. File: {str(img)}")

    def config(self, channels=None, markers=None):
        selected_channels = filter_channels(self, channels=channels)
        config(self, channels=selected_channels, markers=markers)
        self.__stacks = np.asarray([imread(str(self.__channels_files[c])) for c in self.channels])
        return self

    def exp_matrix(self, method="mean", polygonize="convex", alpha=0):
        if len(self.markers) == 0:
            print("ATTENTION: NO marker specific, using channels' name instead.")
        cells = mask2cells(self.__mask_img)
        data = get_cell_exp_stack(self.__stacks, cells, method=method)
        print(f"Detected {len(data)} cells.")

        return data, geom_cells(cells, method=polygonize, alpha=alpha)
