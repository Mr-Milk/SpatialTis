from ast import literal_eval
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import chi2, norm
from tqdm import tqdm

from spatialtis.abc import AnalysisBase
from spatialtis.config import CONFIG
from spatialtis.spatial.utils import QuadStats
from spatialtis.typing import Number, Tuple
from spatialtis.utils import create_remote, doc, run_ray


def VMR(points, resample, r, pval):
    tree = cKDTree(points)
    tmax = tree.maxes
    tmin = tree.mins

    counts = []
    for i in range(0, resample):
        # select a random point
        x = np.random.randint(tmin[0], tmax[0] + 1, 1)
        y = np.random.randint(tmin[1], tmax[1] + 1, 1)
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

        n = len(points)
        chi2_value = (n - 1) * ID
        p_value = 1 - chi2.cdf(chi2_value, n - 1)
        accept_null = p_value > pval

        if accept_null:
            pattern = 1  # random
        elif (not accept_null) & (ID > 1):
            pattern = 3  # clustered
        else:
            pattern = 2  # regular

        return pattern
    else:
        return 0


def QUAD(points, quad, pval):
    counts = QuadStats(points, nx=quad[0], ny=quad[1]).grid_counts()
    quad_count = list(counts.keys())
    # index of dispersion
    n = len(points)
    sum_x = np.sum(quad_count)
    sum_x_sqr = np.sum(np.square(quad_count))
    ID = n * (sum_x_sqr - sum_x) / (sum_x ** 2 - sum_x)
    chi2_value = ID * (sum_x - 1) + n - sum_x
    p_value = 1 - chi2.cdf(chi2_value, n - 1)
    accept_null = p_value > pval

    if accept_null:
        pattern = 1  # random
    elif (not accept_null) & (ID > 1):
        pattern = 3  # clustered
    else:
        pattern = 2  # regular

    return pattern


def NNS(points, pval):
    tree = cKDTree(points)
    tmax = tree.maxes
    tmin = tree.mins

    area = (tmax[0] - tmin[0]) * (tmax[1] - tmin[1])
    r = [tree.query(c, k=[2])[0][0] for c in points]
    n = len(points)
    sum_r = np.sum(r)
    r_A = sum_r / n
    rho = n / area
    r_E = 1 / (2 * np.sqrt(rho))

    # aggregation index R
    R = r_A / r_E
    z_score = (r_A - r_E) / (0.26136 / np.sqrt(n * rho))

    p_value = norm.sf(abs(z_score)) * 2
    accept_null = p_value > pval

    if accept_null:
        pattern = 1  # random
    elif (not accept_null) & (R > 1):
        pattern = 3  # clustered
    else:
        pattern = 2  # regular

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
    | Clark and Evans aggregation index: R | R = 1  |  R < 1  |  R > 1  |
    +--------------------------------------+--------+---------+---------+

    Args:
        data: {adata}
        method: "vmr", "quad", and "nns" (Default: "nns")
        pval: {pval}
        r: Only use when method="vmr", diameter of sample window
        resample: Only use when method="vmr", the number of random permutations to perform
        quad: Only use when method="quad", how to perform rectangle tessellation
        **kwargs: {analysis_kwargs}


    "quad" is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid.

    """

    def __init__(
        self,
        data: AnnData,
        method: str = "nns",
        pval: float = 0.01,
        r: Optional[Number] = 10,
        resample: int = 500,
        quad: Tuple[int, int] = (10, 10),
        **kwargs,
    ):

        if method == "vmr":
            self.method = "Variance-to-mean ratio"
            self._dist_func = VMR
            self._args = [resample, r, pval]
        elif method == "quad":
            self.method = "Quadratic statistic"
            self._dist_func = QUAD
            self._args = [quad, pval]
        else:
            self.method = "Nearest neighbors search"
            self._dist_func = NNS
            self._args = [pval]

        super().__init__(data, task_name="spatial_distribution", **kwargs)

        df = data.obs[self.exp_obs + [self.centroid_key, self.cell_type_key]]
        groups = df.groupby(self.exp_obs)

        patterns = []
        name_tags = []
        type_tags = []

        need_eval = self.is_col_str(self.centroid_key)

        if self.mp:

            dist_mp = create_remote(self._dist_func)

            jobs = []
            for name, group in groups:
                if isinstance(name, str):
                    name = [name]
                for t, tg in group.groupby(self.cell_type_key):
                    if len(tg) > 1:
                        if need_eval:
                            cells = [literal_eval(c) for c in tg[self.centroid_key]]
                        else:
                            cells = [c for c in tg[self.centroid_key]]
                        jobs.append(dist_mp.remote(cells, *self._args))
                        type_tags.append(t)
                        name_tags.append(name)

            patterns = run_ray(jobs, desc="Finding distribution pattern")

        else:
            for name, group in tqdm(
                groups, **CONFIG.pbar(desc="Finding distribution pattern"),
            ):
                if isinstance(name, str):
                    name = [name]
                for t, tg in group.groupby(self.cell_type_key):
                    if len(tg) > 1:
                        if need_eval:
                            cells = [literal_eval(c) for c in tg[self.centroid_key]]
                        else:
                            cells = [c for c in tg[self.centroid_key]]
                        patterns.append(self._dist_func(cells, *self._args))
                        type_tags.append(t)
                        name_tags.append(name)

        dist_patterns = []
        for n, t, p in zip(name_tags, type_tags, patterns):
            dist_patterns.append([*n, t, p])

        self.result = pd.DataFrame(
            data=dist_patterns, columns=self.exp_obs + ["type", "pattern"]
        )
