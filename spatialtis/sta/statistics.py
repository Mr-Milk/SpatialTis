import numpy as np
import pandas as pd
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

    groups = np.unique(sdata[group_by].to_numpy(dtype=str), axis=0)

    types = np.unique(sdata[type_col].to_numpy(dtype=str))

    if selected_type is not None:
        types = [t for t in selected_type if t in types]

    multicol = pd.MultiIndex.from_tuple(groups, names=group_by)

    all_results = pd.DataFrame(columns=multicol, index=types)
    all_results.index.name = 'type'

    for label in all_results.columns:
        query = ' & '.join([f'{x} == "{y}"' for (x, y) in zip(group_by, label)])
        result = sdata.query(query)[type_col]
        c = Counter(result)
        all_results[label] = [c[t] for t in types]

    return all_results


def cell_components(

        data: AnnData,
        type_col: str,
        group_by: list,
        count_base: str,
        selected_type: list = None,
        metrics: str = 'percentage',
        export_key: str = 'cell_components'

):
    all_results = type_counter(data, type_col, group_by, selected_type=selected_type)

    container = {'parameters': {'group_by': group_by,
                                'counted_column': type_col,
                                'types': all_results.index,
                                },
                 'data': all_results}

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
