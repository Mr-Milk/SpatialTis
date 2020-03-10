from collections import OrderedDict
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import norm
from shapely.geometry import MultiPoint

from spatialtis.config import CONFIG
from ..utils import filter_adata


def _hotspot(cells, grid_size, level, pval):
    bbox = MultiPoint(cells).bounds
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    nx = int(width // grid_size)
    ny = int(height // grid_size)

    N = nx * ny

    # the grid must bigger than 3 * 3 so that it can have neighbors
    if N < 9:
        return ['cold' for i in range(0, len(cells))]
    else:
        w_x = width / nx
        h_y = height / ny
        # modify from PySAL pointpats
        dict_id_count = OrderedDict()
        for i in range(ny):
            for j in range(nx):
                dict_id_count[j + i * nx] = 0
        cells_grid_id = []
        for point in cells:
            index_x = int((point[0] - bbox[0]) // w_x)
            index_y = int((point[1] - bbox[1]) // h_y)
            if index_x == nx:
                index_x -= 1
            if index_y == ny:
                index_y -= 1
            id = index_y * nx + index_x
            cells_grid_id.append(id)
            dict_id_count[id] += 1

        quad_count = np.asarray(list(dict_id_count.values())).reshape(nx, ny)
        idx_points = [(i, j) for i in range(0, nx) for j in range(0, ny)]
        hotsquad = []

        tree = cKDTree(idx_points)
        # parameter in equation
        mean_C = np.mean(quad_count)
        sum_c = np.sum(np.square(quad_count.ravel()))
        S = np.sqrt(sum_c / N - mean_C ** 2)
        # There will be some situation when S == 0
        if S == 0:
            return ['cold' for i in range(0, len(cells))]
        else:
            for p in idx_points:
                # neighbors = tree.query_ball_point(p, r=np.sqrt(2))
                neighbors = tree.query_ball_point(p, r=level * np.sqrt(2))
                pp = [idx_points[i] for i in neighbors]
                ix = np.asarray([p[0] for p in pp])
                iy = np.asarray([p[1] for p in pp])
                sum_wc = np.sum(quad_count[ix.min():ix.max(), iy.min():iy.max()].ravel())

                sum_w = len(neighbors)

                U = np.sqrt((N * sum_w - sum_w ** 2) / (N - 1))
                # U == 0 means the neighbor cells is the entire regions
                # meaning the regions is too small so no significant hotspot
                if U == 0:
                    hotsquad.append(False)
                else:
                    # z score for region
                    z = sum_wc - (mean_C * sum_w / (S * U))
                    p_value = norm.sf(np.abs(z))
                    hot = p_value < pval
                    hotsquad.append(hot)

            marker_hot = []
            # print(cells_grid_id)
            for i in cells_grid_id:
                if hotsquad[i]:
                    marker_hot.append('hot')
                else:
                    marker_hot.append('cold')

        return marker_hot


if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`"
        )


    @ray.remote
    def _hotspot_mp(cells, grid_size, level, pval):
        return _hotspot(cells, grid_size, level, pval)


def hotspot(
        adata: AnnData,
        groupby: Union[Sequence, str, None] = None,
        type_col: Optional[str] = None,
        centroid_col: str = 'centroid',
        selected_types: Optional[Sequence] = None,
        search_level: int = 1,
        grid_size: int = 50,
        pval: float = 0.01,
        export_key: str = 'hotspot',
        mp: bool = False,
):
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL

    df = filter_adata(adata, groupby, type_col, centroid_col, selected_types=selected_types, reset_index=False)
    groups = df.groupby(groupby)

    if mp & (CONFIG.OS in ["Linux", "Darwin"]):
        results = []
        indexs = []
        hotcells = []
        for name, group in groups:
            for t, tg in group.groupby(type_col):
                if len(tg) > 1:
                    cells = [eval(c) for c in tg[centroid_col]]
                    results.append(_hotspot_mp.remote(cells, grid_size, search_level, pval))
                    indexs.append(tg.index)
                else:
                    hotcells.append(pd.Series(['cold' for i in range(0, len(tg))], index=tg.index))

        results = ray.get(results)
        for i, hots in enumerate(results):
            hotcells.append(pd.Series(hots, index=indexs[i]))

    else:
        hotcells = []
        for name, group in groups:
            for t, tg in group.groupby(type_col):
                if len(tg) > 1:
                    cells = [eval(c) for c in tg[centroid_col]]
                    hots = _hotspot(cells, grid_size, search_level, pval)
                    hotcells.append(pd.Series(hots, index=tg.index))
                else:
                    hotcells.append(pd.Series(['cold' for i in range(0, len(tg))], index=tg.index))

    hotcells = pd.concat(hotcells)
    adata.obs[export_key] = hotcells
    # Cell map will leave blank if fill with None value
    adata.obs[export_key].fillna("other")
