from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import chi2, norm

from spatialtis.config import CONFIG
from ..utils import df2adata_uns, filter_adata


def _index_of_dispersion(groups, types, type_col, centroid_col, resample, r):
    names = [n for n, _ in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        for t, tg in group.groupby(type_col):
            if len(tg) > 1:
                cells = [eval(c) for c in tg[centroid_col]]
                tree = cKDTree(cells)
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
                if np.mean(counts) != 0:
                    ID = np.var(counts) / np.mean(counts)

                    if ID == 1:
                        pattern = 1  # random
                    elif ID > 1:
                        pattern = 3  # clustered
                    else:
                        pattern = 2  # regular

                    patterns[name][t] = pattern

    return patterns


def _morisita(groups, types, type_col, centroid_col, quad, pval):

    try:
        from pointpats import PointPattern
        from pointpats.quadrat_statistics import QStatistic
    except ImportError:
        raise ImportError("pointpats not installed, try `pip install pointpats`")

    names = [n for n, _ in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        for t, tg in group.groupby(type_col):
            if len(tg) > 1:
                cells = [eval(c) for c in tg[centroid_col]]

                pp = PointPattern(cells)
                counts = QStatistic(
                    pp, shape="rectangle", nx=quad[0], ny=quad[1]
                ).mr.point_location_sta()
                quad_count = list(counts.keys())
                # index of dispersion
                n = len(cells)
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

                patterns[name][t] = pattern

    return patterns


def _ce_aggindex(groups, types, type_col, centroid_col, pval):
    names = [n for n, _ in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        for t, tg in group.groupby(type_col):
            if len(tg) > 1:
                cells = [eval(c) for c in tg[centroid_col]]
                tree = cKDTree(cells)
                tmax = tree.maxes
                tmin = tree.mins

                area = (tmax[0] - tmin[0]) * (tmax[1] - tmin[1])
                r = [tree.query(c, k=[2])[0][0] for c in cells]
                n = len(cells)
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

                patterns[name][t] = pattern

    return patterns


def spatial_distribution(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    centroid_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    method: str = "ID",
    pval: float = 0.01,
    r: Optional[float] = 10,
    resample: int = 50,
    quad: Sequence[int] = (10, 10),
    export: bool = True,
    export_key: str = "spatial_distribution",
    return_df: bool = False,
    overwrite: bool = False,
):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)

     - random (1)
     - regular (2)
     - cluster (3)

    Three methods are provided

     - Index of Dispersion (ID)
     - Morisita’s index of dispersion (MID)
     - Clark and Evans aggregation index (CE)

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
        selected_types: selected cell types you want to count
        method: Options are 'ID', 'MID', or 'CE'
        pval: if smaller than pval, reject null hypothesis (random distribute)
        r: only use when method='ID', diameter of sample window
        resample: only use when method='ID', the number of random permutations to perform
        quad: only use when method='MID', how to perform rectangle tessellation
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field (Default: "spatial_distribution")
        return_df: whether to return an pandas.DataFrame
        overwrite: whether to overwrite if the key existed

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
        adata, groupby, type_col, centroid_col, selected_types=selected_types
    )
    types = pd.unique(df[type_col])
    groups = df.groupby(groupby)

    if method == "ID":
        patterns = _index_of_dispersion(
            groups, types, type_col, centroid_col, resample, r
        )
    elif method == "MID":
        patterns = _morisita(groups, types, type_col, centroid_col, quad, pval)
    elif method == "CE":
        patterns = _ce_aggindex(groups, types, type_col, centroid_col, pval)
    else:
        raise ValueError(
            f"'{method}' No such method, available options are 'ID','MID','CE'."
        )

    results = pd.DataFrame(patterns)
    results = results.rename_axis(index=["Cell type"], columns=groupby).T

    if export:
        df2adata_uns(results, adata, export_key, overwrite)

    if return_df:
        return results
