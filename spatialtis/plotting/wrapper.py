from anndata import AnnData
import pandas as pd
import pickle
from typing import Sequence

from ._bar_plot import stacked_bar


def cell_components(data: AnnData,
                    group_by: Sequence[str],
                    read_key='cell_components',
                    **kwargs
                    ):
    keys = data.uns.keys()
    if read_key not in keys:
        raise (KeyError, "please specific the key using argument 'read_key'")

    container = data.uns[read_key]['data']

    df = pd.DataFrame(eval(container['df']))
    df.index.set_names(container['iname'], inplace=True)
    stacked_bar(df, group_by, **kwargs)

    '''
    
    mindex = pd.MultiIndex.from_tuples(eval(container['index']))
    df = pd.DataFrame(container['m'], index=mindex, columns=container['columns'])
    df.index.names = container['iname']
    stacked_bar(df, group_by, **kwargs)
    '''