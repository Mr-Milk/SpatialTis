import pickle
from typing import Sequence
import numpy as np
import pandas as pd
from anndata import AnnData

from itertools import combinations
from scipy.stats import chisquare

from ._bar_plot import stacked_bar
from ._violin_plot import violin_plot
from ..utils import adata_uns2df
from ._heatmap_sns import heatmap


def cell_components(
        adata: AnnData,
        groupby: Sequence[str],
        key: str = "cell_components",
        **kwargs
):
    df = adata_uns2df(adata, key)
    p = stacked_bar(df, groupby, **kwargs)

    return p


def cell_density(
        adata: AnnData,
        groupby: Sequence[str],
        key: str = "cell_density",
        **kwargs
):
    df = adata_uns2df(adata, key)
    df = pd.DataFrame(df.stack(), columns=["density"])
    groupby = groupby + ["type"]
    p = violin_plot(df, groupby, "density", **kwargs)

    return p


def cell_co_occurrence(
        adata: AnnData,
        groupby: Sequence[str],
        key: str = "cell_co_occurrence",
        **kwargs
):
    df = adata_uns2df(adata, key)
    tdf = df.astype(int).groupby(level=groupby)
    X = []
    for n, g in tdf:
        c = (g.sum() / len(g)).fillna(0)
        data = []
        for comb in combinations(c, 2):
            if (comb[0] == 0) | (comb[1] == 0):
                p = 0
            else:
                p = chisquare(comb).pvalue
            data.append(p)
        X.append(data)
    pdf = (pd.DataFrame(X) > 0.95).astype(int)
    pdf.index = pd.MultiIndex.from_frame(df.index.to_frame(index=False)[groupby].drop_duplicates())
    pdf.columns = pd.MultiIndex.from_arrays(np.asarray([i for i in combinations(df.columns, 2)]).T,
                                            names=['Cell type1', 'Cell type2'])
    p = heatmap(pdf, row_colors=groupby, col_colors=['Cell type1', 'Cell type2'], colorbar_type='categorical',
                categorical_colorbar_text=['Co-occur', 'Non co-occur'], col_cluster=None, row_cluster=None)
    return p
