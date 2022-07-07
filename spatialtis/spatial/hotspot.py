from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import getis_ord, points_bbox
from typing import Tuple, List

from spatialtis.abc import AnalysisBase
from spatialtis.utils import col2adata, doc


@doc
def hotspot(data: AnnData,
            selected_types: List[str] | np.ndarray = None,
            search_level: int = 3,
            quad: Tuple[int, int] = None,
            rect_side: Tuple[float, float] = None,
            pval: float = 0.01,
            export_key: str = None,
            **kwargs,
            ):
    """`Getis-ord hotspot detection <../about/implementation.html#hotspot-detection>`_

    Used to identify cells that cluster together.

    Parameters
    ----------
    data : {adata}
    selected_types : {selected_types}
    search_level : int, default: 3
        How deep the search level to reach.
    quad : tuple of int, default: (10, 10)
        A tuple (X, Y), Use a grid that is X * Y to tessellation your ROI.
    rect_side : tuple of float
        A tuple (X, Y), Use many rectangles with X * Y side to tessellation your ROI.
    pval : {pval}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """

    ab = AnalysisBase(data, display_name="Hotspot analysis", export_key="hotspot", **kwargs)
    ab.check_cell_type()
    if selected_types is not None:
        ab.export_key = f"hotspot_{'_'.join(selected_types)}"
    else:
        selected_types = ab.cell_types
    if export_key is not None:
        ab.export_key = export_key
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
