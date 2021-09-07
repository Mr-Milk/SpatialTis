from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import norm
from spatialtis_core import getis_ord, points_bbox

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import QuadStats
from spatialtis.typing import Array
from spatialtis.utils import col2adata_obs, doc


def _hotspot(cells, bbox, grid_size, level, pval):
    q = QuadStats(cells, bbox, grid_size=grid_size)

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
                    quad_count[ix.min(): ix.max(), iy.min(): iy.max()].ravel()
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
                    p_value = norm.sf(np.abs(z)) * 2
                    hot = p_value < pval
                    hotsquad.append(hot)

            marker_hot = []

            for i in q.cells_grid_id:
                if hotsquad[i]:
                    marker_hot.append("hot")
                else:
                    marker_hot.append("cold")

        return marker_hot


@doc
class hotspot(AnalysisBase):
    """`Getis-ord hotspot detection <../about/implementation.html#hotspot-detection>`_

    Used to identify cells that cluster together.

    Args:
        data: {adata}
        selected_types: {selected_types}
        search_level: How deep the search level to reach
        grid_size: Length of the side of square grid
        pval: {pval}
        kwargs: {analysis_kwargs}

    """

    def __init__(
            self,
            data: AnnData,
            selected_types: Optional[Array] = None,
            search_level: int = 3,
            quad: Optional[Tuple[int, int]] = None,
            rect_side: Optional[Tuple[float, float]] = None,
            pval: float = 0.01,
            **kwargs
    ):
        super().__init__(data, **kwargs)

        hotcells = []
        for roi_name, roi_data in self.roi_iter(desc="Hotspot analysis"):
            points = read_points(roi_data, self.centroid_key)
            bbox = points_bbox(points)
            for t, g in roi_data.groupby(self.cell_type_key):
                cells = read_points(g, self.centroid_key)
                if t in selected_types:
                    hots = getis_ord(cells, bbox, search_level=search_level, quad=quad, rect_side=rect_side)
                    hotcells.append(pd.Series(hots, index=g.index))

        result = pd.concat(hotcells)
        self.data.obs[self.export_key] = result
        # Cell map will leave blank if fill with None value
        self.data.obs[self.export_key].fillna("other", inplace=True)
        arr = self.data.obs[self.export_key].astype("category")
        arr = arr.cat.rename_categories({True: "hot", False: "cold", "other": "other"})
        # Call this to invoke the print
        col2adata_obs(self.data.obs[self.export_key], self.data, self.export_key)
        self.stop_timer()


    @property
    def result(self):
        return self.data.obs[self.exp_obs + [self.cell_type_key, self.export_key]]
