from collections import Counter
from itertools import combinations_with_replacement, product
from typing import List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from milkviz import anno_clustermap, dot_heatmap
from milkviz.utils import mask_triu
from scipy.stats import pearsonr

from spatialtis import Config, get_result
from spatialtis.utils import doc, log_print
from .utils import pairs_to_adj


def count_size_side(pdata, type_order):
    dot_size, dot_hue = {}, {}
    for comb, df in pdata.groupby(['type1', 'type2']):
        count = {1: 0, 0: 0, -1: 0, **Counter(df['relationship'])}
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
        groupby: Optional[List] = None,
        key: str = "cell_interaction",
        type_order: Optional[List[str]] = None,
        order: bool = True,
        plot_value: str = "relationship",
        **plot_options,
):
    pdata = get_result(data, key)
    uni_types = pd.unique(pdata[['type1', 'type2']].to_numpy().flatten())
    if type_order is None:
        type_order = natsorted(uni_types)
    if use == "heatmap":
        # the index of [1::] is to remove the index columns
        pdata = pdata.pivot_table(columns=['type1', 'type2'],
                                  values=plot_value,
                                  index=pdata.index.names[1::],
                                  fill_value=0)
        pdata = pdata[[tuple(c) for c in product(type_order, repeat=2)]]
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
        dot_size, dot_hue = count_size_side(pdata, type_order)
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
        dot_size_data = dot_size.to_numpy(dtype=int)
        dot_hue_data = dot_hue.to_numpy()
        matrix_hue = matrix.to_numpy() if matrix is not None else None
        xticklabels = dot_size.columns
        yticklabels = dot_size.index
        if not order:
            dot_size_data = mask_triu(dot_size_data)
            dot_hue_data = mask_triu(dot_hue_data)
            matrix_hue = mask_triu(matrix_hue) if matrix is not None else None
            xticklabels = xticklabels[::-1]
        options = dict(
            legend_title="Sign' ROI",
            hue_cbar_title="Interaction",
            matrix_cbar_title="Pearson-R",
            hue_cbar_ticklabels=[" - ", "+"],
            matrix_cbar_ticklabels=["-1", " 1"],
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
        type_order: Optional[List[str]] = None,
        order: bool = True,
        **plot_options,
):
    rdata = get_result(data, key)
    dot_size, dot_hue = count_size_side(rdata, type_order)
    dot_size_data = dot_size.to_numpy(dtype=int)
    dot_hue_data = dot_hue.to_numpy()
    xticklabels = dot_size.columns
    yticklabels = dot_size.index

    if not order:
        dot_size_data = mask_triu(dot_size_data)
        dot_hue_data = mask_triu(dot_hue_data)
        xticklabels = xticklabels[::-1]
    options = dict(
        size_legend_title="Sign' ROI",
        hue_cbar_title="Interaction",
        hue_cbar_ticklabels=[" - ", "+"],
        sizes=(0, 500),
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
