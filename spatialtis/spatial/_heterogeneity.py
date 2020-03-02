from typing import Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from scipy.stats import entropy

from ..sta.statistics import type_counter
from ..utils import df2adata_uns
from spatialtis.config import CONFIG


def spatial_heterogeneity(
        adata: AnnData,
        groupby: Union[Sequence, str, None] = None,
        type_col: Optional[str] = None,
        compare: Optional[int] = None,
        selected_types: Optional[Sequence] = None,
        export_key: str = 'spatial_heterogeneity',
        return_df: bool = False,
) -> Optional[pd.DataFrame]:
    """compute spatial heterogeneity
    Here we use entropy for spatial heterogeneity, which describes the amount of information.
    To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences
    for each sample within the group are computed, smaller value indicates less difference within group.

    Args:
        adata: anndata object to perform analysis
        groupby: list of names describes your experiments
        type_col: name of the cell types column
        compare: Compute Kullback-Leibler divergences based on which level
        selected_types: which cell type to select, I highly recommended you not to set this args, just let it be None.
        export_key: the key name to store info, exported to anndata.uns field
        return_df: whether to return a pandas DataFrame object

    """
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL

    df = type_counter(adata, groupby, type_col, selected_types=selected_types)

    if len(df.columns) == 1:
        print("No heterogeneity, only one type of cell found.")
        return None

    KL_div = {}
    if compare is not None:
        groups = df.groupby(level=groupby[compare])
        for n, g in groups:
            count_g = g.sum()
            qk = count_g.div(count_g.sum())
            KL_div[n] = qk

    ent = []
    KL = []
    KL_level = []
    for row in df.iterrows():
        pk = list(row[1].div(row[1].sum()))
        if compare is not None:
            compare_level = row[0][compare]
            KL.append(entropy(pk, KL_div[compare_level], base=2))
            KL_level.append(compare_level)
        ent.append(entropy(pk, base=2))

    data = {'heterogeneity': ent}
    if compare is not None:
        data['KL'] = KL
        data['level'] = KL_level
    roi_heterogeneity = pd.DataFrame(data=data, index=df.index)

    # export to anndata
    df2adata_uns(roi_heterogeneity, adata, export_key, overwrite='True')

    if return_df:
        return roi_heterogeneity
