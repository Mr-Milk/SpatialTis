from collections import OrderedDict
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import norm
from tqdm import tqdm

from spatialtis.config import CONFIG

from ..utils import filter_adata
from ._util import quad_sta


def _hotspot(cells, grid_size, level, pval):
    q = quad_sta(cells, grid_size=grid_size)

    nx = q.nx
    ny = q.ny
    N = nx * ny

    # the grid must bigger than 3 * 3 so that it can have neighbors
    if N < 9:
        return ["cold" for _ in range(0, len(cells))]
    else:
        dict_id_count = q.grid_counts()

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
            return ["cold" for _ in range(0, len(cells))]
        else:
            for p in idx_points:
                # neighbors = tree.query_ball_point(p, r=np.sqrt(2))
                neighbors = tree.query_ball_point(p, r=level * np.sqrt(2))
                pp = [idx_points[i] for i in neighbors]
                ix = np.asarray([p[0] for p in pp])
                iy = np.asarray([p[1] for p in pp])
                sum_wc = np.sum(
                    quad_count[ix.min() : ix.max(), iy.min() : iy.max()].ravel()
                )

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

            for i in q.cells_grid_id:
                if hotsquad[i]:
                    marker_hot.append("hot")
                else:
                    marker_hot.append("cold")

        return marker_hot


if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`",
        )

    _hotspot_mp = ray.remote(_hotspot)


def hotspot(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    search_level: int = 1,
    grid_size: int = 50,
    pval: float = 0.01,
    export_key: Optional[str] = None,
    mp: bool = False,
):
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    if export_key is None:
        export_key = CONFIG.hotspot_key
    else:
        CONFIG.hotspot_key = export_key

    df = filter_adata(
        adata,
        groupby,
        type_key,
        centroid_key,
        selected_types=selected_types,
        reset_index=False,
    )
    groups = df.groupby(groupby)

    if mp & (CONFIG.OS in ["Linux", "Darwin"]):

        def exec_iterator(obj_ids):
            while obj_ids:
                done, obj_ids = ray.wait(obj_ids)
                yield ray.get(done[0])

        results = []
        indexs = []
        hotcells = []
        for name, group in groups:
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [eval(c) for c in tg[centroid_key]]
                    results.append(
                        _hotspot_mp.remote(cells, grid_size, search_level, pval)
                    )
                    indexs.append(tg.index)
                elif len(tg) == 1:
                    hotcells.append(pd.Series(["cold"], index=tg.index))

        for _ in tqdm(
            exec_iterator(results),
            total=len(results),
            desc="hotspot analysis",
            bar_format=CONFIG.PBAR_FORMAT,
            disable=(not CONFIG.PROGRESS_BAR),
        ):
            pass

        results = ray.get(results)

        for i, hots in enumerate(results):
            hotcells.append(pd.Series(hots, index=indexs[i]))

    else:
        hotcells = []
        for name, group in tqdm(
            groups, bar_format=CONFIG.PBAR_FORMAT, disable=(not CONFIG.PROGRESS_BAR)
        ):
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [eval(c) for c in tg[centroid_key]]
                    hots = _hotspot(cells, grid_size, search_level, pval)
                    hotcells.append(pd.Series(hots, index=tg.index))
                elif len(tg) == 1:
                    hotcells.append(pd.Series(["cold"], index=tg.index))

    hotcells = pd.concat(hotcells)
    adata.obs[export_key] = hotcells
    # Cell map will leave blank if fill with None value
    adata.obs[export_key].fillna("other")
