from anndata import AnnData

from typing import Sequence, Union, Optional

from ..sta.statistics import type_counter


def spatial_heterogeneity(
        adata: AnnData,
        groupby: Union[Sequence, str],
        type_col: str,
        selected_types: Optional[Sequence] = None,

):

    df = type_counter(adata, groupby, type_col, selected_types=selected_types)

