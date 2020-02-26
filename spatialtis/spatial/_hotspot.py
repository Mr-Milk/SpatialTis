import numpy as np
import pandas as pd
from anndata import AnnData

from shapely.geometry import MultiPoint
from scipy.spatial import cKDTree
from scipy.stats import norm
from pointpats import PointPattern
from pointpats.quadrat_statistics import QStatistic

from typing import Sequence, Union, Optional
from collections import OrderedDict

from ..utils import filter_adata
from ._util import quad_tessellation


def hotspot(
        adata: AnnData,
        groupby: Union[Sequence, str],
        type_col: str,
        centroid_col: str = 'centroid',
        selected_types: Optional[Sequence] = None,
        search_level: int = 3,
        grid_size: int = 50,
        pval: float = 0.01,
        export_key: str = 'hotspot'
):
    df = filter_adata(adata, groupby, type_col, centroid_col, selected_types=selected_types, reset_index=False)
    groups = df.groupby(groupby)
    hotcells = []
    for name, group in groups:
        for t, tg in group.groupby(type_col):
            cells = tg[centroid_col]
            hotcells += _hotspot(cells, grid_size, search_level, pval)

    df['hotspot'] = hotcells
    adata.obs[export_key] = df['hotspot']


def _hotspot(cells, grid_size, level, pval):
    bbox = MultiPoint(cells).bounds
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    nx = width // grid_size
    ny = height // grid_size

    if nx < 1:
        nx = 1
    if ny < 1:
        ny = 1

    # modify from PySAL pointpats
    dict_id_count = OrderedDict()
    for i in range(ny):
        for j in range(nx):
            dict_id_count[j + i * nx] = 0
    cells_grid_id = []
    for point in cells:
        index_x = (point[0] - bbox[0]) // width
        index_y = (point[1] - bbox[1]) // height
        if index_x == nx:
            index_x -= 1
        if index_y == ny:
            index_y -= 1
        id = index_y * nx + index_x
        cells_grid_id.append(id)
        dict_id_count[id] += 1

    quad_count = dict_id_count.values()
    idx_points = [(i, j) for i in range(0, nx) for j in range(0, ny)]
    hotsquad = []

    tree = cKDTree(idx_points)
    # parameter in equation
    N = nx * ny
    mean_C = np.mean(quad_count)
    for p in idx_points:
        neighbors = tree.query_ball_point(p, r=np.sqrt(2))
        all_points = tree.query_ball_point(p, r=level * np.sqrt(2))
        pp = [idx_points[i] for i in neighbors]
        ix = np.asarray([p[0] for p in pp])
        iy = np.asarray([p[1] for p in pp])
        sum_wc = np.sum(quad_count[ix.min():ix.max(), iy.min():iy.max()])

        pp = [idx_points[i] for i in all_points]
        ix = np.asarray([p[0] for p in pp])
        iy = np.asarray([p[1] for p in pp])
        sum_c = np.sum(quad_count[ix.min():ix.max(), iy.min():iy.max()])

        sum_w = len(neighbors)

        # z score for region
        S = np.sqrt(sum_c / N - mean_C ** 2)
        U = np.sqrt((N * sum_w - sum_w ** 2) / (N - 1))
        z = sum_wc - (mean_C * sum_w / (S * U))
        p_value = 1 - norm.cdf(z)
        accept_null = p_value > pval
        hotsquad.append(accept_null)

    marker_hot = []
    for i in cells_grid_id:
        if hotsquad[i]:
            marker_hot.append(None)
        else:
            marker_hot.append('hot')

    return marker_hot
