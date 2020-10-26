from ast import literal_eval
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import chi2, norm
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.utils import QuadStats
from spatialtis.utils import (
    create_remote,
    df2adata_uns,
    filter_adata,
    get_default_params,
    log_print,
    reuse_docstring,
    run_ray,
    timer,
)


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


@timer(prefix="Running spatial distribution")
@get_default_params
@reuse_docstring()
def spatial_distribution(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    method: str = "nns",
    pval: float = 0.01,
    r: Optional[float] = 10,
    resample: int = 500,
    quad: Sequence[int] = (10, 10),
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)

     - random (1)
     - regular (2)
     - cluster (3)

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
        adata: {adata}
        groupby: {groupby}
        method: Options are "vmr", "quad", and "nns"
        pval: {pval}
        r: Only use when method="vmr", diameter of sample window
        resample: Only use when method="vmr", the number of random permutations to perform
        quad: Only use when method="quad", how to perform rectangle tessellation
        type_key: {type_key}
        centroid_key: {centroid_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}

    "quad" is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid.

    """
    if export_key is None:
        export_key = CONFIG.spatial_distribution_key
    else:
        CONFIG.spatial_distribution_key = export_key

    df = filter_adata(adata, groupby, type_key, centroid_key,)
    types = pd.unique(df[type_key])
    groups = df.groupby(groupby)

    if method == "vmr":
        log_print("Method: Variance-to-mean ratio")
        _dist_func = VMR
        args = [resample, r, pval]
    elif method == "quad":
        log_print("Method: Quadratic statistic")
        _dist_func = QUAD
        args = [quad, pval]
    elif method == "nns":
        log_print("Method: Nearest neighbors search")
        _dist_func = NNS
        args = [pval]
    else:
        raise ValueError(
            f"'{method}' No such method, available options are 'vmr','quad' and 'nns'."
        )

    patterns = []
    name_tags = []
    type_tags = []

    if mp:

        _dist_mp = create_remote(_dist_func)

        jobs = []
        for name, group in groups:
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [literal_eval(c) for c in tg[centroid_key]]
                    jobs.append(_dist_mp.remote(cells, *args))
                    type_tags.append(t)
                    name_tags.append(name)

        patterns = run_ray(
            jobs,
            tqdm_config=CONFIG.tqdm(
                total=len(jobs), desc="Finding distribution pattern", unit="task"
            ),
        )

    else:
        for name, group in tqdm(
            groups, **CONFIG.tqdm(desc="Finding distribution pattern"),
        ):
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [literal_eval(c) for c in tg[centroid_key]]
                    patterns.append(_dist_func(cells, *args))
                    type_tags.append(t)
                    name_tags.append(name)

    dist_patterns = {n: {t: 0 for t in types} for n in name_tags}
    for n, t, p in zip(name_tags, type_tags, patterns):
        dist_patterns[n][t] = p

    results = pd.DataFrame(dist_patterns)
    results = results.rename_axis(index=["Cell type"], columns=groupby).T

    if export:
        df2adata_uns(results, adata, export_key, params={"method": method})

    if return_df:
        return results
