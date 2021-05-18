from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import chi2, norm
from shapely.geometry import MultiPoint

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import QuadStats, get_eval
from spatialtis.typing import Number, Tuple
from spatialtis.utils import create_remote, doc, run_ray
from spatialtis.utils.log import pbar_iter


def get_pattern(ID, pvalue, pval):
    reject_null = pvalue < pval

    if reject_null:
        if ID > 1:
            pattern = 3  # cluster
        elif ID == 1:
            pattern = 2  # regular
        else:
            pattern = 1  # random
    else:
        pattern = 1  # random
    return pattern


def VMR(points, bbox, min_cells, pval, resample, r):
    n = len(points)
    if n < min_cells:
        return 0
    else:
        minx, miny, maxx, maxy = bbox
        tree = cKDTree(points)

        counts = []
        for i in range(0, resample):
            # select a random point
            x = np.random.randint(minx, maxx + 1, 1)
            y = np.random.randint(miny, maxy + 1, 1)
            # query the point
            query_point = [x[0], y[0]]
            neighbor_points = tree.query_ball_point(query_point, r)
            counts.append(len(neighbor_points))

        # index of dispersion
        counts = np.array(counts)
        # since this is a sampling method, there is still a very small probability
        # that we sample nothing
        if np.mean(counts) != 0:
            ID = np.var(counts) / np.mean(counts)
            chi2_value = (n - 1) * ID
            p_value = 1 - chi2.cdf(chi2_value, n - 1)
            pattern = get_pattern(ID, p_value, pval)
            return pattern
        else:
            return 0


def QUAD(points, bbox, min_cells, pval, quad=None, grid_size=None):
    n = len(points)
    if n < min_cells:
        return 0
    else:
        if quad is not None:
            counts = QuadStats(points, bbox, nx=quad[0], ny=quad[1]).grid_counts()
        else:
            counts = QuadStats(points, bbox, grid_size=grid_size).grid_counts()
        quad_count = np.asarray(list(counts.keys()))
        # index of dispersion
        sum_x = np.sum(quad_count)
        sum_x_sqr = np.sum(np.square(quad_count))
        if sum_x > 1:
            ID = n * (sum_x_sqr - sum_x) / (sum_x ** 2 - sum_x)
            chi2_value = ID * (sum_x - 1) + n - sum_x
            p_value = 1 - chi2.cdf(chi2_value, n - 1)
            pattern = get_pattern(ID, p_value, pval)
        else:
            # when there is only one cell or no cells in the grid
            # it will cause ZeroDivision error
            pattern = 0
        return pattern


def NNS(points, bbox, min_cells, pval):
    n = len(points)
    if n < min_cells:
        return 0
    else:
        minx, miny, maxx, maxy = bbox
        tree = cKDTree(points)

        area = (maxx - minx) * (maxy - miny)
        r = np.array([tree.query(c, k=[2])[0][0] for c in points])
        intensity = n / area
        # sum_r = np.sum(r)
        # r_A = sum_r / n
        nnd_mean = r.mean()
        nnd_expected_mean = 1 / (2 * np.sqrt(intensity))
        R = nnd_mean / nnd_expected_mean

        SE = np.sqrt(((4 - np.pi) * area) / (4 * np.pi)) / n
        Z = (nnd_mean - nnd_expected_mean) / SE

        p_value = norm.sf(abs(Z)) * 2
        reject_null = p_value < pval

        if reject_null:
            if R < 1:
                pattern = 3
            elif R == 1:
                pattern = 2
            else:
                pattern = 1
        else:
            pattern = 1  # random
        return pattern


@doc
class spatial_distribution(AnalysisBase):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)

    - Random (1)
    - Regular (2)
    - Cluster (3)

    Three methods are provided

    - Variance-to-mean ratio (vmr): `Index of Dispersion <../about/implementation.html#index-of-dispersion>`_
    - Quadratic statistics (quad): `Morisita’s index of dispersion <../about/implementation.html#morisitas-index-of-dispersion>`_
    - Nearest neighbors search (nns): `Clark and Evans aggregation index <../about/implementation.html#clark-and-evans-aggregation-index>`_

    +--------------------------------------+--------+---------+---------+
    |                                      | Random | Regular | Cluster |
    +======================================+========+=========+=========+
    | Index of dispersion: ID              | ID = 1 | ID < 1  | ID > 1  |
    +--------------------------------------+--------+---------+---------+
    | Morisita’s index of dispersion: I    | I = 1  |  I < 1  |  I > 1  |
    +--------------------------------------+--------+---------+---------+
    | Clark and Evans aggregation index: R | R = 1  |  R > 1  |  R < 1  |
    +--------------------------------------+--------+---------+---------+

    Args:
        data: {adata}
        method: "vmr", "quad", and "nns" (Default: "nns")
        min_cells: The minimum number of the specific type of cells in a ROI to perform analysis
        pval: {pval}
        r: Only use when method="vmr", determine diameter of sample window, should be in [0, 1], default is 0.1
            this take 1/10 of the shortest side of the ROI as the diameter.
        resample: Only use when method="vmr", the number of random permutations to perform
        quad: Only use when method="quad", how to perform rectangle tessellation. Default is (10, 10), this will use a
            10*10 grid to perform tessellation.
        grid_size: Only use when method="quad", the side of grid when perform rectangle tessellation.
        **kwargs: {analysis_kwargs}


    "quad" is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid.

    """

    def __init__(
            self,
            data: AnnData,
            method: str = "nns",
            min_cells: int = 5,
            pval: float = 0.01,
            r: Number = 0.1,
            resample: int = 500,
            quad: Optional[Tuple[int, int]] = None,
            grid_size: Optional[Number] = None,
            **kwargs,
    ):

        if method == "vmr":
            self.method = "Variance-to-mean ratio"
            self._dist_func = VMR
        elif method == "quad":
            self.method = "Quadratic statistic"
            self._dist_func = QUAD
            if quad is not None:
                self._args = [quad]
            else:
                if grid_size is not None:
                    self._args = [None, grid_size]
                else:
                    self._args = [(10, 10)]
        else:
            self.method = "Nearest neighbors search"
            self._dist_func = NNS
            self._args = []

        super().__init__(data, task_name="spatial_distribution", **kwargs)

        df = data.obs[self.exp_obs + [self.centroid_key, self.cell_type_key]]
        groups = df.groupby(self.exp_obs)
        need_eval = self.is_col_str(self.centroid_key)

        patterns = []
        name_tags = []
        type_tags = []

        if self.mp:

            dist_mp = create_remote(self._dist_func)

            jobs = []
            for name, group in groups:
                if isinstance(name, str):
                    name = [name]
                ROI = get_eval(group, self.centroid_key, need_eval)
                bbox = MultiPoint(ROI).bounds

                if method == "vmr":
                    # auto generate a r parameters for every ROI
                    # 1/10 of the shortest side
                    minx, miny, maxx, maxy = bbox
                    roi_r = min([maxx - minx, maxy - miny]) * r
                    self._args = [resample, roi_r]

                for t, tg in group.groupby(self.cell_type_key, sort=False):
                    cells = get_eval(tg, self.centroid_key, need_eval)
                    jobs.append(
                        dist_mp.remote(cells, bbox, min_cells, pval, *self._args)
                    )
                    name_tags.append(name)
                    type_tags.append(t)

            patterns = run_ray(jobs, desc="Finding distribution pattern")

        else:
            for name, group in pbar_iter(
                    groups, desc="Finding distribution pattern",
            ):
                if isinstance(name, str):
                    name = [name]
                ROI = get_eval(group, self.centroid_key, need_eval)
                bbox = MultiPoint(ROI).bounds

                if method == "vmr":
                    minx, miny, maxx, maxy = bbox
                    roi_r = min([maxx - minx, maxy - miny]) * r
                    self._args = [resample, roi_r]

                for t, tg in group.groupby(self.cell_type_key, sort=False):
                    cells = get_eval(tg, self.centroid_key, need_eval)
                    patterns.append(
                        self._dist_func(cells, bbox, min_cells, pval, *self._args)
                    )
                    name_tags.append(name)
                    type_tags.append(t)

        dist_patterns = []
        for n, t, p in zip(name_tags, type_tags, patterns):
            dist_patterns.append([*n, t, p])
        self.result = pd.DataFrame(
            data=dist_patterns, columns=self.exp_obs + ["type", "pattern"]
        )
