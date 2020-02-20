import numpy as np
import pandas as pd
from anndata import AnnData

from typing import Sequence, Union, Optional

from ..utils import filter_adata


def hotspot(
        adata: AnnData,
        groupby: Union[Sequence, str],
        type_col: str,
        centroid_col: str = 'centroid',
        selected_types: Optional[Sequence] = None,
        pval: float = 0.01,
):
    df = filter_adata(adata, groupby, type_col, centroid_col, selected_types=selected_types, reset_index=False)
    groups = df.groupby(groupby)
    for name, group in groups:
        pass

