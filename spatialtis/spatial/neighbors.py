import warnings
from ast import literal_eval
from typing import Optional

import pandas as pd
from anndata import AnnData

from spatialtis.abc import AnalysisBase
from spatialtis.config import CONFIG
from spatialtis.typing import Number
from spatialtis.utils import col2adata_obs, doc
from spatialtis.utils.log import pbar_iter


@doc
class find_neighbors(AnalysisBase):
    """To `find the neighbors <../about/implementation.html#find-cell-neighbors>`_ of each cell

    KD-tree is used when cells are points

    R-tree is used when cells are polygons (use_shape=True)

    Args:
        data: {adata}
        expand: If the cell has shape, it means how much units to expand each cell;
                If the cell is point, it's the search radius;
        scale: How much to scale each cell, only if cell has shape
        count: Count the number of neighbors for each cell (Default: True)
        use_shape: Cell representation is shape instead of point (Default: False)
        **kwargs: {analysis_kwargs}

    """

    def __init__(
            self,
            data: AnnData,
            expand: Optional[Number] = None,
            scale: Optional[Number] = None,
            count: Optional[bool] = True,
            use_shape: Optional[bool] = False,
            **kwargs,
    ):
        if use_shape:
            self.method = "R-tree"
        else:
            self.method = "KD-tree"
        if "export_key" not in kwargs:
            kwargs["export_key"] = CONFIG.NEIGHBORS_KEY
        super().__init__(data, task_name="find_neighbors", **kwargs)

        try:
            import neighborhood_analysis as na
        except ImportError:
            raise ImportError(
                "Package not found, try `pip install neighborhood_analysis`."
            )

        if (expand is None) & (scale is None):
            raise ValueError("Neither 'expand' or 'scale' are specific")

        if expand is not None:
            if expand < 0:
                raise ValueError("Can't shrink cell, 'expand' must >= 0")

        if scale is not None:
            if scale < 1:
                raise ValueError("Can't shrink cell, 'scale' must >= 1")

        if (expand is not None) & (scale is not None):
            warnings.warn(
                f"Conflict parameters, can't set 'expand' and 'scale' in the same time, use expand={expand}"
            )

        if (expand is None) & (not use_shape):
            raise ValueError("Cannot scale a point, use 'expand' instead")

        # prepare data for shape neighbor search
        self.data.obs[self.neighbors_ix_key] = [
            i for i in range(self.data.obs.shape[0])
        ]
        track_ix = []
        neighbors = []
        if use_shape:
            need_eval = self.is_col_str(self.shape_key)
            bboxs = []
            labels = []
            for n, g in pbar_iter(
                    data.obs.groupby(self.exp_obs, sort=False),
                    desc="Get cell bbox",
            ):
                shapes = g[self.shape_key]
                if need_eval:
                    cell_shapes = [literal_eval(cell) for cell in shapes]
                else:
                    cell_shapes = [cell for cell in shapes]
                bbox = na.get_bbox(cell_shapes)
                bboxs.append(bbox)
                labels.append(list(g[self.neighbors_ix_key]))
                track_ix += list(g.index)

            for bbox, label in pbar_iter(
                    zip(bboxs, labels), desc="Find neighbors", total=len(bbox),
            ):
                if expand is not None:
                    neighbors += na.get_bbox_neighbors(
                        bbox, expand=expand, labels=label
                    )
                else:
                    neighbors += na.get_bbox_neighbors(bbox, scale=scale, labels=label)
        else:
            groups = data.obs.groupby(self.exp_obs, sort=False)
            for n, g in pbar_iter(
                    groups,
                    desc="Find neighbors",
                    total=len(groups),
            ):
                cells = [literal_eval(c) for c in g[self.centroid_key]]
                neighbors += na.get_point_neighbors(
                    cells, expand, labels=list(g[self.neighbors_ix_key])
                )
                track_ix += list(g.index)
        col2adata_obs(pd.Series(neighbors, index=track_ix), data, self.export_key)
        if count:
            counts = [len(i) for i in neighbors]
            col2adata_obs(counts, data, f"{self.export_key}_count")
        self.stop_timer()

    @property
    def result(self):
        return self.data.obs[self.exp_obs + [self.export_key]]
