import warnings
from pathlib import Path
from typing import Optional, Sequence, Union

import anndata as ad
import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.preprocessing.geom import get_cell_exp_stack, mask2cells
from spatialtis.utils import create_remote, reuse_docstring, run_ray


def get_roi(
    exp_img,
    mask_img,
    bg: Optional[int] = 0,
    method: str = "mean",
    polygonize: str = "convex",
    alpha: Optional[float] = None,
):
    # From skimage doc: The different color bands/channels are stored in the third dimension
    # so we need to transpose it
    # exp_img = np.transpose(imread(str(exp_img)), (2, 1, 0))

    # read page by page, we don't know how user will store their file
    # some store 40 channels on one page, some store 1 channel per page
    try:
        import tifffile
    except ImportError:
        raise ImportError("Required tifffile, try `pip install tifffile`.")

    exp = []
    with tifffile.TiffFile(str(exp_img)) as img:
        for i in img.pages:
            exp.append(i.asarray())

    exp = np.asarray(exp)
    if exp.ndim > 4:
        raise ValueError("The dimensions of image are too much")
    else:
        # if the channels info store in the
        if exp.shape[0] == 1:
            exp = np.transpose(exp[0], (2, 1, 0))

    cells, geom_info = mask2cells(mask_img, bg=bg, polygonize=polygonize, alpha=alpha)
    data = get_cell_exp_stack(exp, cells, method=method)

    return [data, geom_info]


class read_ROIs:
    """Extract single cell expression matrix and geometry information from stacked images and masks

    Args:
        entry: The root folder to start with
        obs_names: Array of names correspond to each level of your folders
        var: Describe the order of layers in your stacked image
        mask_pattern: Name pattern for all of your mask
        img_pattern: Name pattern for all of your image

    Attributes:
        obs: Will pass to `AnnData.obs`
        var: Will pass to `AnnData.var`
        anndata: The processed `AnnData` object

    """

    def __init__(
        self,
        entry: Union[Path, str],
        obs_names: Sequence,
        var: pd.DataFrame,
        mask_pattern: Optional[str] = None,
        img_pattern: Optional[str] = None,
    ):
        self._obs_names = obs_names
        self._tree = []
        self._mask_img = []
        self._exp_img = []
        self._exhaust_dir(entry)

        self.anndata = None
        self.obs = []
        self.var = var
        # anndata require str index, hard set everything to str
        self.var.index = self.var.index.map(str)

        # make sure every end-dir are at the same level
        level = None
        obs_count = len(obs_names)
        for t in self._tree:
            parts = t.parts
            if level is None:
                level = len(parts)

            if level == len(parts):
                self.obs.append(parts[-obs_count:])
                # locate the mask image
                mask_set = [mask for mask in t.glob(f"*{mask_pattern}*")]
                if len(mask_set) > 1:
                    raise ValueError(f"More than one mask image found, {t}")
                elif len(mask_set) == 0:
                    raise ValueError(f"No mask image found, {t}")
                else:
                    self._mask_img.append(mask_set[0])

                # locate the exp image
                exp_set = [exp for exp in t.glob(f"*{img_pattern}*")]
                if len(exp_set) > 1:
                    raise ValueError(
                        f"More than one image found, "
                        "please stacked them together and delete the unused. {t}"
                    )
                elif len(exp_set) == 0:
                    raise ValueError(f"No image found, {t}")
                else:
                    self._exp_img.append(exp_set[0])
            else:
                raise ValueError("The depth of your file directory are not consistent")

    # walk through the directory, until there is no directory
    def _exhaust_dir(
        self, path: Union[Path, str],
    ):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self._tree.append(f)
            if f.parent in self._tree:
                self._tree.remove(f.parent)
            self._exhaust_dir(f)

    @reuse_docstring()
    def to_anndata(
        self,
        bg: Optional[int] = 0,
        method: str = "mean",
        polygonize: str = "convex",
        alpha: Optional[float] = None,
        mp: Optional[bool] = None,
    ):
        """Get anndata object

        You must explicitly call this method to trigger the computation.

        Args:
            bg: The background pixel value
            method: How to compute the expression level. ("mean", "sum", "median")
            polygonize: How to compute the cell shape.("convex", "concave")
            alpha: The alpha value for polygonize="concave"
            mp: {mp}

        .. note:: **"convex" or "concave" to determine cell shape?**

                    The cell shape is represent by the border points to simplify the following analysis process.

                    - **convex**: Convex hull, much faster but less accurate.
                    - **concave**: Concave hull, very slow, a parameter "alpha" is needed.


        """

        if mp is None:
            mp = CONFIG.MULTI_PROCESSING

        X = []
        ann_obs = []

        areas = []
        shapes = []
        centroids = []
        eccentricities = []

        if polygonize == "concave":
            warnings.warn("Running concave hull is very slow", RuntimeWarning)

        if mp:

            get_roi_mp = create_remote(get_roi)

            jobs = []
            for exp_img, mask_img in zip(self._exp_img, self._mask_img):
                jobs.append(
                    get_roi_mp.remote(
                        exp_img,
                        mask_img,
                        bg=bg,
                        method=method,
                        polygonize=polygonize,
                        alpha=alpha,
                    )
                )

            mp_results = run_ray(
                jobs,
                tqdm_config=CONFIG.tqdm(total=len(self._tree), desc="process images"),
            )

            for (exp, cells), obs in zip(mp_results, self.obs):
                X += exp
                ann_obs += list(np.repeat(np.array([obs]), len(cells[0]), axis=0))
                areas += cells[0]
                shapes += cells[1]
                centroids += cells[2]
                eccentricities += cells[3]

        else:
            for exp_img, mask_img, obs in tqdm(
                zip(self._exp_img, self._mask_img, self.obs),
                **CONFIG.tqdm(total=len(self._tree), desc="process images"),
            ):
                [exp, cells] = get_roi(
                    exp_img,
                    mask_img,
                    bg=bg,
                    method=method,
                    polygonize=polygonize,
                    alpha=alpha,
                )
                X += exp
                ann_obs += list(np.repeat(np.array([obs]), len(cells[0]), axis=0))
                areas += cells[0]
                shapes += cells[1]
                centroids += cells[2]
                eccentricities += cells[3]

        # print(len(ann_obs), len(areas))
        # anndata require str index, hard set to str
        ann_obs = pd.DataFrame(
            ann_obs,
            columns=self._obs_names,
            index=[str(i) for i in range(0, len(ann_obs))],
        )
        ann_obs[CONFIG.AREA_KEY] = areas
        ann_obs[CONFIG.SHAPE_KEY] = shapes
        ann_obs[CONFIG.CENTROID_KEY] = centroids
        ann_obs[CONFIG.ECCENTRICITY_KEY] = eccentricities

        X = np.asarray(X, dtype=float)

        self.anndata = ad.AnnData(X, obs=ann_obs, var=self.var, dtype="float")

        return self.anndata
