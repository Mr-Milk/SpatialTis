import numpy as np
import pandas as pd
import pickle
from collections import Counter

from anndata import AnnData


def type_counter(
        data: AnnData,
        type_col: str,
        group_by: list = None,
        selected_type: list = None
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

    return all_results  # obs as index, var as columns


def cell_components(

        data: AnnData,
        type_col: str,
        group_by: list,
        selected_type: list = None,
        export_key: str = 'cell_components'

):
    all_results = type_counter(data, type_col, group_by, selected_type=selected_type)

    container = {'parameters': {'group_by': group_by,
                                'counted_column': type_col,
                                'types': list(all_results.columns),
                                },
                 'data': dict(df=str(all_results.to_dict()),
                              iname=all_results.index.names,
                              columns=list(all_results.columns))}

    # write to anndata object
    data.uns[export_key] = container
    print(f'''Finished!
    Add to AnnData object
    uns: '{export_key}' ''')


def cell_co_occurrence(

        data: AnnData,
        type_col: str,
        group_by: list = None,
        selected_type: list = None,
        threshold: int = 50,
        export_key: str = 'co_occurrence'

):
    all_results = type_counter(data, type_col, group_by, selected_type=selected_type)

    occurrence = pd.DataFrame(index=all_results.index)
    for col, con in all_results.items():
        occur = con > threshold
        occurrence[col] = occur

    data.uns[export_key] = occurrence
    print(f'''Finished!
    Add to AnnData object
    uns: '{export_key}' ''')
