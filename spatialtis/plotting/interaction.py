from collections import Counter

import numpy as np
import pandas as pd
from anndata import AnnData
from itertools import product
from milkviz import anno_clustermap, dot_heatmap
from milkviz.utils import mask_triu
from natsort import natsorted
from scipy.stats import pearsonr
from typing import List

from spatialtis import get_result
from spatialtis.utils import doc, log_print, df2adata_uns
from .utils import pairs_to_adj


def count_size_side(pdata, type_order, groupby_keys, value_key):
    dot_size, dot_hue = {}, {}
    for comb, df in pdata.groupby(groupby_keys):
        count = {1: 0, 0: 0, -1: 0, **Counter(df[value_key])}
        v = [1, -1]
        arr = [count[i] for i in v]
        sig_count = np.sum(arr)
        if sig_count == 0:
            sig_num = 0
        else:
            norm = arr / sig_count
            sig_num = v[np.argmax(norm)] * np.amax(norm)
        dot_size[comb] = sig_count
        dot_hue[comb] = sig_num
    dot_size = pairs_to_adj(
        pd.DataFrame(dot_size, index=[0]).T.reset_index(), type_order
    )
    dot_hue = pairs_to_adj(pd.DataFrame(dot_hue, index=[0]).T.reset_index(), type_order)
    return dot_size, dot_hue


def count_size_side_for_enrichment(pdata, type_order):
    dot_size, dot_hue = {}, {}
    for comb, arr in pdata.iteritems():
        count = {1: 0, 0: 0, -1: 0, **Counter(arr)}
        v = [1, -1]
        arr = [count[i] for i in v]
        sig_count = np.sum(arr)
        if sig_count == 0:
            sig_num = 0
        else:
            norm = arr / sig_count
            sig_num = v[np.argmax(norm)] * np.amax(norm)
        dot_size[comb] = sig_count
        dot_hue[comb] = sig_num
    dot_size = pairs_to_adj(
        pd.DataFrame(dot_size, index=[0]).T.reset_index(), type_order
    )
    dot_hue = pairs_to_adj(pd.DataFrame(dot_hue, index=[0]).T.reset_index(), type_order)
    return dot_size, dot_hue


@doc
def cell_interaction(
        data: AnnData,
        use: str = "dot",
        groupby: List = None,
        key: str = "cell_interaction",
        type_order: List[str] = None,
        order: bool = True,
        plot_value: str = "relationship",
        **plot_options,
):
    """Visualization of the cell interaction analysis

    Parameters
    ----------
    data : {adata_plotting}
    use : {'dot', 'heatmap'}, default: 'dot'
    groupby : {groupby}
    key : {plot_key}
    type_order: {type_order}
    order : bool
    plot_value : {'relationship', 'statistic'}
    **plot_options :
        Pass to :func:`milkviz.dot_heatmap` or :func:`milkviz.anno_clustermap`.

    """
    if use == "heatmap":

        store_key = "interaction_heatmap"
        if store_key in data.uns_keys():
            pdata = get_result(data, store_key)
        else:
            pdata = get_result(data, key)
            uni_types = pd.unique(pdata[['type1', 'type2']].to_numpy().flatten())
            if type_order is None:
                type_order = natsorted(uni_types)
            pdata = pdata.pivot_table(columns=['type1', 'type2'],
                                      values=plot_value,
                                      # the index of [1::] is to remove the index columns
                                      index=pdata.index.names[1::],
                                      fill_value=0)
            pdata = pdata[[tuple(c) for c in product(type_order, repeat=2)]]
            df2adata_uns(pdata, data, store_key, verbose=False)

        if plot_value == "relationship":
            options = dict(
                categorical_cbar=["Avoidance", "Association"],
                col_legend_split=False,
                col_legend_title="Cell type",
                cbar_title="Interaction",
                col_cluster=False,
                method="ward",
                vmin=-1,
                vmax=1,
            )
            options = {**options, **plot_options}
            return anno_clustermap(
                pdata, col_colors=["type1", "type2"], row_colors=groupby, **options
            )
        else:
            options = dict(
                col_legend_split=False,
                col_legend_title="Cell type",
                cbar_title="Interaction",
                col_cluster=False,
                method="ward",
            )
            options = {**options, **plot_options}
            return anno_clustermap(
                pdata, col_colors=["type1", "type2"], row_colors=groupby, **options
            )
    else:
        def _get_cell_components():
            try:
                matrix = get_result(data, "cell_components")
                combs = {}
                for i1, i2 in product(matrix.iteritems(), repeat=2):
                    combs[(i1[0], i2[0])] = pearsonr(i1[1], i2[1])[0]
                matrix = pd.DataFrame(combs, index=[0]).T.reset_index()
                matrix = pairs_to_adj(matrix, type_order)
                matrix = matrix.loc[
                    dot_size.index, dot_size.columns
                ]  # ensure the number match to data
            except Exception:
                log_print(
                    "Run spatialtis.cell_components to add "
                    "pearson's R of cell proportion to the visualization"
                )
                matrix = None
            return matrix

        incomplete = False
        save = False
        store_key = "interaction_dot"
        if store_key in data.uns_keys():
            pdata = data.uns[store_key]
            dot_size_data = pdata['dot_size_data']
            dot_hue_data = pdata['dot_hue_data']
            matrix_hue = pdata['matrix_hue']
            xticklabels = pdata['xticklabels']
            yticklabels = pdata['yticklabels']
            if pdata['incomplete']:
                incomplete = False
                save = True
                matrix_hue = _get_cell_components()
        else:
            save = True
            pdata = get_result(data, key)
            uni_types = pd.unique(pdata[['type1', 'type2']].to_numpy().flatten())
            if type_order is None:
                type_order = natsorted(uni_types)
            dot_size, dot_hue = count_size_side(pdata, type_order,
                                                groupby_keys=['type1', 'type2'], value_key="relationship")
            matrix = _get_cell_components()
            if matrix is None:
                matrix_hue = None
                incomplete = True
            else:
                matrix_hue = matrix.to_numpy()
            dot_size_data = dot_size.to_numpy(dtype=int)
            dot_hue_data = dot_hue.to_numpy()
            xticklabels = dot_size.columns
            yticklabels = dot_size.index
            if not order:
                dot_size_data = mask_triu(dot_size_data)
                dot_hue_data = mask_triu(dot_hue_data)
                matrix_hue = mask_triu(matrix_hue) if matrix is not None else None
                xticklabels = xticklabels[::-1]

        if save:
            pdata = dict(
                dot_size_data=dot_size_data,
                dot_hue_data=dot_hue_data,
                matrix_hue=matrix_hue,
                xticklabels=xticklabels,
                yticklabels=yticklabels,
                incomplete=incomplete
            )
            data.uns[store_key] = pdata

        options = dict(
            dot_size_legend_kw={"title": "Sign' ROI"},
            dot_hue_cbar_kw={"title": "Interaction"},
            matrix_cbar_kw={"title": "Pearson-R"},
            sizes=(0, 250),
            dot_cmap="RdYlBu_r",
        )
        options = {**options, **plot_options}
        return dot_heatmap(
            dot_size=dot_size_data,
            dot_hue=dot_hue_data,
            matrix_hue=matrix_hue,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            **options,
        )


@doc
def spatial_enrichment(
        data: AnnData,
        key: str = "spatial_enrichment",
        type_order: List[str] = None,
        **plot_options,
):
    """Visualization of the spatial enrichment analysis

    Parameters
    ----------
    data : {adata_plotting}
    key : {plot_key}
    type_order : {type_order}
    **plot_options :
        Pass to :func:`milkviz.dot_heatmap`.

    """
    store_key = "spatial_enrichment_dot"
    if store_key in data.uns_keys():
        pdata = data.uns[store_key]
        dot_size_data = pdata['dot_size_data']
        dot_hue_data = pdata['dot_hue_data']
        xticklabels = pdata['xticklabels']
        yticklabels = pdata['yticklabels']
    else:
        rdata = get_result(data, key)
        dot_size, dot_hue = count_size_side_for_enrichment(rdata, type_order)
        dot_size_data = dot_size.to_numpy(dtype=int)
        dot_hue_data = dot_hue.to_numpy()
        xticklabels = dot_size.columns
        yticklabels = dot_size.index

        pdata = dict(
            dot_size_data=dot_size_data,
            dot_hue_data=dot_hue_data,
            xticklabels=xticklabels,
            yticklabels=yticklabels
        )
        data.uns[store_key] = pdata

    # if not order:
    #     dot_size_data = mask_triu(dot_size_data)
    #     dot_hue_data = mask_triu(dot_hue_data)
    #     xticklabels = xticklabels[::-1]
    options = dict(
        dot_size_legend_kw={"title": "Sign' ROI"},
        dot_hue_cbar_kw={"title": "Interaction"},
        matrix_cbar_kw={"title": "Pearson-R"},
        sizes=(0, 250),
        dot_cmap="PiYG_r",
        **plot_options,
    )
    return dot_heatmap(
        dot_size=dot_size_data,
        dot_hue=dot_hue_data,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        **options,
    )
