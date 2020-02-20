import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import chisquare, chi2, norm
from anndata import AnnData
from shapely.geometry import box, Point

from typing import Sequence, Union, Optional

from ..utils import filter_adata


def _index_of_dispersion(groups, types, type_col, centroid_col, resample, r, pval):
    names = [n for n in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        typ = group[type_col]
        ct = group[centroid_col]

        uq_typ = pd.unique(typ)
        ctdb = {t: [] for t in uq_typ}

        for i, point in enumerate(ct):
            ctdb[typ[i]].append(point)

        for t, points in ctdb.items():
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
            # p-value
            p = chisquare(counts, [1 for i in range(0, len(counts))]).pvalue

            # index of dispersion
            counts = np.array(counts)
            ID = np.var(counts) / np.mean(counts)
            accept_null = p > pval

            if accept_null:
                pattern = 1  # random
            elif (not accept_null) & (ID > 1):
                pattern = 3  # clustered
            else:
                pattern = 2  # regular

            patterns[name][t] = pattern

    return patterns


def _rec_tesllation(
        rect: tuple,
        quad: Union[tuple, str] = 'auto'
):
    min_unit = 150
    if quad == 'auto':
        nx = int((rect[2] - rect[0]) / 20)
        ny = int((rect[3] - rect[1]) / 20)
        quad = (nx, ny)

    xrange = np.linspace(rect[0], rect[2], quad[0] + 1)
    yrange = np.linspace(rect[1], rect[3], quad[1] + 1)

    rects = []
    for i, x in enumerate(xrange[0:-1]):
        x1 = x
        x2 = xrange[i + 1]
        lower = [(x1, y) for y in yrange[0:-1]]
        higher = [(x2, y) for y in yrange[1::]]
        for p in tuple(zip(lower, higher)):
            rects.append(box(*tuple(np.asarray(p).ravel())))
    return rects


def _morisita(groups, types, type_col, centroid_col, quad, pval):
    names = [n for n in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        typ = group[type_col]
        ct = group[centroid_col]

        uq_typ = pd.unique(typ)
        ctdb = {t: [] for t in uq_typ}

        x_coord = [p[0] for p in ct]
        y_coord = [p[1] for p in ct]
        xmin = np.min(x_coord)
        ymin = np.min(y_coord)
        xmax = np.max(x_coord)
        ymax = np.max(y_coord)

        rects = _rec_tesllation((xmin, ymin, xmax, ymax), quad=quad)

        for i, point in enumerate(ct):
            ctdb[typ[i]].append(point)

        for t, points in ctdb.items():
            rect_quads = [0 for i in range(0, len(rects))]
            for p in points:
                for i, rect in enumerate(rects):
                    if rect.contains(Point(p)):
                        rect_quads[i] += 1

            # index of dispersion
            n = len(rects)
            sum_x = np.sum(rect_quads)
            sum_x_sqr = np.sum(np.square(rect_quads))
            ID = n * (sum_x_sqr - sum_x) / (sum_x ** 2 - sum_x)
            chi2_value = ID * (sum_x - 1) + n - sum_x
            p_value = 1 - chi2.cdf(chi2_value, n - 1)
            accept_null = p_value < pval

            if accept_null:
                pattern = 1  # random
            elif (not accept_null) & (ID > 1):
                pattern = 3  # clustered
            else:
                pattern = 2  # regular

            patterns[name][t] = pattern

    return patterns


def _ce_aggindex(groups, types, type_col, centroid_col, pval):
    names = [n for n in groups]
    patterns = {n: {t: 0 for t in types} for n in names}

    for name, group in groups:
        typ = group[type_col]
        ct = group[centroid_col]

        uq_typ = pd.unique(typ)
        ctdb = {t: [] for t in uq_typ}

        for i, point in enumerate(ct):
            ctdb[typ[i]].append(point)

        for t, points in ctdb.items():
            tree = cKDTree(points)
            tmax = tree.maxes
            tmin = tree.mins

            area = (tmax[0] - tmin[0])*(tmax[1] - tmin[1])
            r = [tree.query(p)[0] for p in points]
            n = len(points)
            sum_r = np.sum(r)
            r_A = sum_r / n
            rho = n / area
            r_E = 1 / (2*np.sqrt(rho))

            # aggregation index R
            R = r_A / r_E
            z_score = (r_A - r_E) / (0.26136 / np.sqrt(n * rho))

            p_value = 1 - norm.cdf(z_score)
            accept_null = p_value < pval

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
        groupby: Union[Sequence, str],
        type_col: str,
        centroid_col: str = 'centroid',
        selected_types: Optional[Sequence] = None,
        method: str = 'ID',
        pval: float = 0.01,
        r: Optional[float] = None,
        resample: int = 50,
        quad: Union[tuple, str] = 'auto',
):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)
    - random (1)
    - regular (2)
    - clumped (3)

    Three methods are provided
    - Index of Dispersion (ID)
    - Morisita’s index of dispersion (MID)
    - Clark and Evans aggregation index (CE)

    Args:
        adata: anndata object to perform analysis
        groupby: list of name describe your experiments
        type_col: name of the cell types column
        centroid_col: name of the cell's centroid column
        selected_types: which cell type
        method: Use 'ID', 'MID', or 'CE'
        pval: if smaller than pval, reject null hypothesis (random distribute)
        r: only use when method='ID', diameter of sample window
        resample: only use when method='ID', the number of random permutations to perform
        quad: only use when method='MID', how to perform rectangle tessellation

    MID is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid,
    if you don't know what to choose, let us choose for you, the default is 'auto'
    """
    df = filter_adata(adata, groupby, type_col, centroid_col, selected_types=selected_types)
    types = df.unique(df[type_col])
    groups = df.groupby(groupby)

    if method == 'ID':
        patterns = _index_of_dispersion(groups, types, type_col, centroid_col, resample, r, pval)
    elif method == 'MID':
        patterns = _morisita(groups, types, type_col, centroid_col, quad, pval)
    elif method == 'CE':
        patterns = _ce_aggindex(groups, types, type_col, centroid_col, pval)
    else:
        raise ValueError(f"'{method}' No such method, available options are 'ID','MID','CE'.")

    results = pd.DataFrame(patterns)

    return results
