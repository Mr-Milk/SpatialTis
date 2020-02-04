from anndata import AnnData
import pandas as pd
import pickle
from typing import Sequence

from ._bar_plot import stacked_bar
from ._violin_plot import violin_plot


def _plot_df(data: AnnData, key=None):
    keys = data.uns.keys()
    if key not in keys:
        raise KeyError(f"{key} not found, please specific using 'key' argument")

    container = data.uns[key]
    df = pd.DataFrame(eval(container['df']))
    df.index.set_names(container['iname'], inplace=True)

    return df


def cell_components(data: AnnData,
                    group_by: Sequence[str],
                    key='cell_components',
                    **kwargs
                    ):
    df = _plot_df(data, key)
    p = stacked_bar(df, group_by, **kwargs)

    return p


def cell_density(data: AnnData,
                 group_by: Sequence[str],
                 key='cell_density',
                 **kwargs
                 ):

    df = _plot_df(data, key)
    df = pd.DataFrame(df.stack(), columns=['density'])
    group_by = group_by + ['type']
    p = violin_plot(df, group_by, 'density', **kwargs)

    return p
