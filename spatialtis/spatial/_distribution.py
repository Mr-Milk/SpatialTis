from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from tqdm import tqdm
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import chi2, norm

from spatialtis.config import CONFIG
from ..utils import df2adata_uns, filter_adata
from ._util import quad_sta


def _wrapper(groups, types, type_col, centroid_col, patch_func, *args):
    names = [n for n, _ in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in tqdm(groups, unit="ROI", desc="find distribution pattern",
                            bar_format=CONFIG.PBAR_FORMAT, disable=(not CONFIG.PROGRESS_BAR)):
        for t, tg in group.groupby(type_col):
            if len(tg) > 1:
                cells = [eval(c) for c in tg[centroid_col]]
                pattern = patch_func(cells, *args)
                patterns[name][t] = pattern

    return patterns


def VMR(points, resample, r):

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

        if ID == 1:
            pattern = 1  # random
        elif ID > 1:
            pattern = 3  # clustered
        else:
            pattern = 2  # regular

        return pattern
    else:
        return 0


def QUAD(points, quad, pval):

    counts = quad_sta(points, nx=quad[0], ny=quad[1]).grid_counts()
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

    p_value = 1 - norm.cdf(z_score)
    accept_null = p_value > pval

    if accept_null:
        pattern = 1  # random
    elif (not accept_null) & (R > 1):
        pattern = 3  # clustered
    else:
        pattern = 2  # regular

    return pattern


def spatial_distribution(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    centroid_col: Optional[str] = None,
    method: str = "nns",
    pval: float = 0.01,
    r: Optional[float] = 10,
    resample: int = 50,
    quad: Sequence[int] = (10, 10),
    export: bool = True,
    export_key: str = "spatial_distribution",
    return_df: bool = False,
):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)

     - random (1)
     - regular (2)
     - cluster (3)

    Three methods are provided

     - Variance-to-mean ratio (vmr): Index of Dispersion (ID)
     - Quadratic statistics (quad): Morisita’s index of dispersion (MID)
     - Nearest neighbors search (nns): Clark and Evans aggregation index (CE)

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
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
        centroid_col: anndata.obs key that store cell centroid info
        method: Options are 'vmr', 'quad', or 'nns'
        pval: if smaller than pval, reject null hypothesis (random distribute)
        r: only use when method='vmr', diameter of sample window
        resample: only use when method='vmr', the number of random permutations to perform
        quad: only use when method='quad', how to perform rectangle tessellation
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field (Default: "spatial_distribution")
        return_df: whether to return an pandas.DataFrame

    MID is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid,
    if you don't know what to choose, let us choose for you, the default is 'auto'

    """
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL

    df = filter_adata(
        adata, groupby, type_col, centroid_col,
    )
    types = pd.unique(df[type_col])
    groups = df.groupby(groupby)

    if method == "vmr":
        patterns = _wrapper(groups, types, type_col, centroid_col, VMR, resample, r)
    elif method == "quad":
        patterns = _wrapper(groups, types, type_col, centroid_col, QUAD, quad, pval)
    elif method == "nns":
        patterns = _wrapper(groups, types, type_col, centroid_col, NNS, pval)
    else:
        raise ValueError(
            f"'{method}' No such method, available options are 'vmr','quad','nns'."
        )

    results = pd.DataFrame(patterns)
    results = results.rename_axis(index=["Cell type"], columns=groupby).T

    if export:
        df2adata_uns(results, adata, export_key)

    if return_df:
        return results
