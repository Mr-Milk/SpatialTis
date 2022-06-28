from typing import Optional, Tuple

import pandas as pd
from anndata import AnnData
from spatialtis_core import getis_ord, points_bbox

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Array
from spatialtis.utils import col2adata, doc


@doc
def hotspot(data: AnnData,
            selected_types: Optional[Array] = None,
            search_level: int = 3,
            quad: Optional[Tuple[int, int]] = None,
            rect_side: Optional[Tuple[float, float]] = None,
            pval: float = 0.01,
            export_key: str = "hotspot",
            **kwargs, ):
    """`Getis-ord hotspot detection <../about/implementation.html#hotspot-detection>`_

    Used to identify cells that cluster together.

    Args:
        data: {adata}
        selected_types: {selected_types}
        search_level: How deep the search level to reach
        quad: {quad}
        rect_side: {rect_size}
        pval: {pval}
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    """

    ab = AnalysisBase(data, display_name="Hotspot analysis", export_key="hotspot_all", **kwargs)
    ab.check_cell_type()
    if selected_types is not None:
        ab.export_key = f"hotspot_{'_'.join(selected_types)}"
    else:
        selected_types = ab.cell_types
    hotcells = []
    for roi_name, cell_types, points, ix in ab.iter_roi(fields=['cell_type', 'centroid', 'index']):
        bbox = points_bbox(points)
        roi_iter = pd.DataFrame({
            "cells": points,
            "types": cell_types
        }, index=ix)
        for t, g in roi_iter.groupby("types"):
            cells = g['cells']
            if t in selected_types:
                hots = getis_ord(
                    cells,
                    bbox,
                    search_level=search_level,
                    quad=quad,
                    rect_side=rect_side,
                    pval=pval,
                )
                hotcells.append(pd.Series(hots, index=g.index))

    result = pd.concat(hotcells)
    data.obs[ab.export_key] = result
    # Cell map will leave blank if fill with None value
    data.obs[ab.export_key].fillna("other", inplace=True)
    arr = data.obs[ab.export_key].astype("category")
    arr = arr.cat.rename_categories({True: "hot", False: "cold", "other": "other"})
    data.obs[ab.export_key] = arr
    # Call this to invoke the print
    col2adata(data.obs[ab.export_key], data, ab.export_key)
    ab.stop_timer()
