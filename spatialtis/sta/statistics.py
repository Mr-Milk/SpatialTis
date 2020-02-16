import numpy as np
import pandas as pd
import pickle
from collections import Counter

from anndata import AnnData

from typing import Sequence, Optional, Union

Num = Union[int, float]


def type_counter(
        data: AnnData,
        type_col: str,
        group_by: Union[Sequence, str],
        selected_type: Optional[Sequence] = None
) -> pd.DataFrame:

    if isinstance(group_by, str):
        sdata = data.obs[[group_by, type_col]]
    else:
        sdata = data.obs[list(group_by) + [type_col]]

    groups = sdata.groupby(group_by)
    types = np.unique(sdata[type_col].to_numpy(dtype=str))

    if selected_type is not None:
        types = [t for t in selected_type if t in types]

    matrix = list()
    mindex = list()
    for i, (n, g) in enumerate(groups):
        c = Counter(g.leiden)
        matrix.append([c[t] for t in types])
        if isinstance(n, str):
            mindex.append((n, i,))
        else:
            mindex.append((*n, i,))

    multiIndex = pd.MultiIndex.from_tuples(mindex)
    all_results = pd.DataFrame(dict(zip(types, np.array(matrix).T)), index=multiIndex)
    all_results.index.set_names(group_by + ['id'], inplace=True)
    all_results.columns.set_names('type', inplace=True)

    return all_results  # obs as index, var as columns


def cell_components(

        data: AnnData,
        type_col: str,
        group_by: list,
        selected_type: Optional[Sequence] = None,
        export: bool = True,
        export_key: str = 'cell_components',
        return_df: bool = False,

):
    results = type_counter(data, type_col, group_by, selected_type=selected_type)

    # Currently, anndata from bioconda is < 0.7.0, pypi is updated
    # It haven't implemented storage of dataframe
    # plus, it mess up with multiIndex dataframe

    if export:
        container = dict(df=str(results.to_dict()),
                         iname=results.index.names,
                         colname=list(results.columns.names))
        # write to anndata object
        data.uns[export_key] = container
        print(f'''Finished!
        Add to AnnData object
        uns: '{export_key}' ''')

    if return_df:
        return results


def cell_co_occurrence(

        data: AnnData,
        type_col: str,
        group_by: list,
        selected_type: list = None,
        export: bool = True,
        export_key: str = 'co_occurrence',
        threshold: int = 50,

):
    counter = type_counter(data, type_col, group_by, selected_type=selected_type)

    occurrence = counter.gt(threshold)

    data.uns[export_key] = occurrence
    print(f'''Finished!
    Add to AnnData object
    uns: '{export_key}' ''')


def cell_density(
        data: AnnData,
        type_col: str,
        group_by: list,
        size: Union[Sequence[Sequence[Num]], Sequence[Num]],
        selected_type: Optional[Sequence] = None,
        export: bool = True,
        export_key: str = 'cell_density',
        ratio: float = 1.0,
        return_df: bool = False,
):
    counter = type_counter(data, type_col, group_by, selected_type=selected_type)
    if isinstance(size[0], (int, float)):
        area = size[0] * size[1]
        results = counter.div(area * ratio)

    elif isinstance(size[0], Sequence):
        area = pd.Series({'area': [s[0] * s[1] * ratio for s in size]})
        results = counter.div(area['area'], axis=0)

    if export:
        container = dict(df=str(results.to_dict()),
                         iname=results.index.names,
                         colname=list(results.columns.names))
        # write to anndata object
        data.uns[export_key] = container
        print(f'''Finished!
        Add to AnnData object
        uns: '{export_key}' ''')

    if return_df:
        return results
