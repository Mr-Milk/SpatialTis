import pickle
from typing import Sequence

import pandas as pd
from anndata import AnnData

from ._bar_plot import stacked_bar
from ._violin_plot import violin_plot
from ..utils import adata_uns2df


def cell_components(
    data: AnnData, groupby: Sequence[str], key="cell_components", **kwargs
):
    df = adata_uns2df(data, key)
    p = stacked_bar(df, groupby, **kwargs)

    return p


def cell_density(data: AnnData, groupby: Sequence[str], key="cell_density", **kwargs):

    df = adata_uns2df(data, key)
    df = pd.DataFrame(df.stack(), columns=["density"])
    groupby = groupby + ["type"]
    p = violin_plot(df, groupby, "density", **kwargs)

    return p
