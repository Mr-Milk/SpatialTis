import warnings
from pathlib import Path
from typing import Optional, Sequence, Union

import anndata as ad
import numpy as np
import pandas as pd

from ._utils import config, config_file, read_ROI, set_info
from spatialtis.config import CONFIG

if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`"
        )


    @ray.remote
    def _get_roi(t, channels, markers, pg, mt, obsi, stacked):
        if len(markers) >= 1:
            roi = read_ROI(t, stacked=stacked).config(
                channels=channels, markers=markers
            )
        else:
            roi = read_ROI(t, stacked=stacked).config(channels=channels)

        exp, cells = roi.exp_matrix(polygonize=pg, method=mt)

        cell_count = len(cells[0])
        obs = np.repeat(np.array([obsi]), cell_count, axis=0)
        print(f"Added: {' '.join(obsi)}")
        return [exp, list(obs), cells]


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

        self.__stacked = (stacked,)
        self.__mask_pattern = mask_pattern
        self.__channels_files = dict()

        self.tree = []
        self.__obs_name = conditions
        self.__depth = len(conditions)

        self.__exhaust_dir(entry)
        self.obs = [p.parts[-self.__depth :] for p in self.tree]

        self.var = None

    # walk through the directory, until there is no directory
    def __exhaust_dir(
        self, path: Union[Path, str],
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
                raise ValueError(
                    "If using concave method, alpha value should > 0, please set argument `alpha`"
                )
            warnings.warn("Running alphashape is very slow", RuntimeWarning)
        elif polygonize == "convex":
            pass
        else:
            raise ValueError("Polygonize options are 'convex' or 'concave'")

        if mp & CONFIG.OS in ["Linux", "Darwin"]:

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
                        self.__stacked,
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
                    roi = read_ROI(d, stacked=self.__stacked).config(
                        channels=self.channels, markers=self.markers
                    )
                else:
                    roi = read_ROI(d, stacked=self.__stacked).config(
                        channels=self.channels
                    )

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