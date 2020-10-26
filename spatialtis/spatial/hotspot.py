from ast import literal_eval
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import norm
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.utils import QuadStats
from spatialtis.utils import (
    col2adata_obs,
    create_remote,
    filter_adata,
    get_default_params,
    reuse_docstring,
    run_ray,
    timer,
)


def _hotspot(cells, grid_size, level, pval):
    q = QuadStats(cells, grid_size=grid_size)

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


@timer(prefix="Running hotspot detection")
@get_default_params
@reuse_docstring()
def hotspot(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    search_level: int = 1,
    grid_size: int = 50,
    pval: float = 0.01,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
):
    """`Getis-ord hotspot detection <../about/implementation.html#hotspot-detection>`_

    Args:
        adata: {adata}
        groupby: {groupby}
        selected_types: {selected_types}
        search_level: How deep the search level to reach
        grid_size: Length of the side of square grid
        pval: {pval}
        type_key: {type_key}
        centroid_key: {centroid_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}

    """
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

    hotcells = []
    if mp:

        _hotspot_mp = create_remote(_hotspot)

        jobs = []
        indexs = []
        for name, group in groups:
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [literal_eval(c) for c in tg[centroid_key]]
                    jobs.append(
                        _hotspot_mp.remote(cells, grid_size, search_level, pval)
                    )
                    indexs.append(tg.index)
                elif len(tg) == 1:
                    hotcells.append(pd.Series(["cold"], index=tg.index))

        results = run_ray(
            jobs,
            tqdm_config=CONFIG.tqdm(
                total=len(jobs), desc="Hotspot analysis", unit="task"
            ),
        )

        for hots, i in zip(results, indexs):
            hotcells.append(pd.Series(hots, index=i))

    else:
        for name, group in tqdm(groups, **CONFIG.tqdm(desc="Hotspot analysis")):
            for t, tg in group.groupby(type_key):
                if len(tg) > 1:
                    cells = [literal_eval(c) for c in tg[centroid_key]]
                    hots = _hotspot(cells, grid_size, search_level, pval)
                    hotcells.append(pd.Series(hots, index=tg.index))
                elif len(tg) == 1:
                    hotcells.append(pd.Series(["cold"], index=tg.index))

    hotcells = pd.concat(hotcells)

    if export:
        adata.obs[export_key] = hotcells
        # Cell map will leave blank if fill with None value
        adata.obs[export_key].fillna("other", inplace=True)
        col2adata_obs(adata.obs[export_key], adata, export_key)

    if return_df:
        return hotcells
