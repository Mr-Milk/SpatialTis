from pathlib import Path
from typing import Optional, Sequence, Union

import anndata as ad
import numpy as np
import pandas as pd
from spatialtis_core import dumps_wkt_points, dumps_wkt_polygons, points_shapes

from spatialtis.config import Config
from spatialtis.utils import create_remote, doc, pbar_iter, run_ray


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
        from tifffile import TiffFile
        from skimage.io import imread
        from skimage.measure import label, regionprops
    except ImportError:
        raise ImportError("Required scikit-image, try `pip install scikit-image`.")

    exp = []
    with TiffFile(str(exp_img)) as img:
        for i in img.pages:
            exp.append(i.asarray())

    exp = np.asarray(exp)
    X, Y, C = 1, 2, 0
    if exp.ndim > 4:
        raise ValueError("The dimensions of image are too much")
    else:
        # if the channels info store in one page
        if exp.shape[0] == 1:
            exp = exp[0]
            X, Y, C = 0, 1, 2
    borders, centroids = [], []
    cells = []
    mask = imread(mask_img)
    label_mask = label(mask, background=bg)
    for cell in regionprops(label_mask):
        # if cell has less than 3 points, then it's meaningless
        if len(cell.coords) >= 3:
            border = points_shapes([(x, y) for x, y in cell.coords], method=polygonize, concavity=alpha)
            borders.append(border)
            centroids.append(cell.centroid)
            cells.append(cell.coords)
    if C == 0:
        cell_exp = [
            getattr(np, method).__call__(exp[:, [i[0] for i in cell], [i[1] for i in cell]], axis=1)
            for cell in cells
        ]
    else:
        cell_exp = [
            getattr(np, method).__call__(exp[[i[0] for i in cell], [i[1] for i in cell], :], axis=0)
            for cell in cells
        ]

    return cell_exp, borders, centroids


class read_ROIs:
    """Extract single cell expression matrix and geometry information from stacked images and masks

    For details, please refer to :ref:`Multiplexed images to AnnData`.

    >>> import spatialtis as st

    >>> var = pd.read_csv("panel_markers.csv")

    >>> reader = st.read_ROIs(entry="IMC_images",
                              obs_names=['ROI'],
                              var=var,
                              mask_pattern="mask",
                              img_pattern="full")
    >>> data = reader.to_anndata()


    Args:
        entry: The root folder to start with
        obs_names: Array of names correspond to each level of your folders
        var: Describe the order of layers in your stacked image
        mask_pattern: Name pattern for all of your mask
        img_pattern: Name pattern for all of your image

    Attributes:
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
        self,
        path: Union[Path, str],
    ):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self._tree.append(f)
            if f.parent in self._tree:
                self._tree.remove(f.parent)
            self._exhaust_dir(f)

    @doc
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

        mp = Config.mp if mp is None else mp

        X, ann_obs, shapes, centroids = [], [], [], []

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

            mp_results = run_ray(jobs, desc="Process images")

            for (exp, borders, centroids_), obs in zip(
                mp_results, self.obs
            ):
                X += exp
                cell_count = len(exp)
                ann_obs += list(np.repeat(np.array([obs]), cell_count, axis=0))
                shapes += borders
                centroids += centroids_

        else:
            for exp_img, mask_img, obs in pbar_iter(
                zip(self._exp_img, self._mask_img, self.obs),
                desc="Process images",
                total=len(self._exp_img),
            ):
                exp, borders, centroids_ = get_roi(
                    exp_img,
                    mask_img,
                    bg=bg,
                    method=method,
                    polygonize=polygonize,
                    alpha=alpha,
                )
                X += exp
                cell_count = len(exp)
                ann_obs += list(np.repeat(np.array([obs]), cell_count, axis=0))
                shapes += borders
                centroids += centroids_

        # anndata require str index, hard set to str
        ann_obs = pd.DataFrame(
            ann_obs,
            columns=self._obs_names,
            index=[str(i) for i in range(0, len(ann_obs))],
        )

        ann_obs["cell_shape"] = dumps_wkt_polygons(shapes)
        ann_obs["centroid"] = dumps_wkt_points(centroids)

        X = np.asarray(X, dtype=float)
        self.anndata = ad.AnnData(X, obs=ann_obs, var=self.var, dtype="float")

        return self.anndata
