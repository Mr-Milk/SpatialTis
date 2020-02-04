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
        group_by: list,
        selected_type: Optional[Sequence] = None
) -> pd.DataFrame:
    sindex = group_by.copy()
    sindex.append(type_col)
    sdata = data.obs[sindex]

    groups = sdata.groupby(group_by)
    types = np.unique(sdata[type_col].to_numpy(dtype=str))

    if selected_type is not None:
        types = [t for t in selected_type if t in types]

    matrix = list()
    mindex = list()
    for n, g in groups:
        c = Counter(g.leiden)
        matrix.append([c[t] for t in types])
        mindex.append(n)

    multiIndex = pd.MultiIndex.from_tuples(mindex)
    all_results = pd.DataFrame(dict(zip(types, np.array(matrix).T)), index=multiIndex)
    all_results.index.set_names(group_by, inplace=True)
    all_results.columns.set_names('type', inplace=True)

    return all_results  # obs as index, var as columns


def cell_components(

        data: AnnData,
        type_col: str,
        group_by: list,
        selected_type: Optional[Sequence] = None,
        export: bool = True,
        export_key: str = 'cell_components',

):
    results = type_counter(data, type_col, group_by, selected_type=selected_type)

    # Currently, anndata from bioconda is < 0.7.0, pypi is updated
    # It haven't implemented storage of dataframe
    # plus, it mess up with multiIndex dataframe

    if export:
        container = dict(df=str(results.to_dict()),
                         iname=results.index.names,
                         columns=list(results.columns))
        # write to anndata object
        data.uns[export_key] = container
        print(f'''Finished!
        Add to AnnData object
        uns: '{export_key}' ''')

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
        ratio: float = 1.0
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
                         columns=list(results.columns))
        # write to anndata object
        data.uns[export_key] = container
        print(f'''Finished!
        Add to AnnData object
        uns: '{export_key}' ''')

    return results
