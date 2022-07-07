from collections import Counter

import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib.colors import ListedColormap
from milkviz import dot, anno_clustermap
from natsort import natsorted
from typing import List

from spatialtis import get_result
from spatialtis.utils import doc


@doc
def spatial_heterogeneity(
        data: AnnData,
        groupby: str = None,
        key: str = "heterogeneity",
        **plot_options,
):
    """Visualization of spatial heterogeneity analysis

    Parameters
    ----------
    data : {adata_plotting}
    groupby : {groupby}
    key : {plot_key}
    **plot_options :
        Pass to :func:`seaborn.boxplot`.

    """
    pdata = get_result(data, key).reset_index()
    exp_obs = pdata.columns[0:-1].tolist()
    if (len(exp_obs) == 1) | (groupby == exp_obs[0]):
        options = dict(color="#5DAC81", **plot_options)
        ax = sns.barplot(data=pdata, y="heterogeneity", x=exp_obs[-1], **options)
    else:
        groupby = exp_obs[0] if groupby is None else groupby
        ax = sns.boxplot(data=pdata, y="heterogeneity", x=groupby, **plot_options)
    return ax


@doc
def cell_dispersion(
        data: AnnData,
        use: str = "dot",
        groupby: List[str] = None,
        type_order: List[str] = None,
        key: str = "cell_dispersion",
        **plot_options,
):
    """

    Parameters
    ----------
    data : {adata_plotting}
    use : {'dot', 'heatmap'}
    groupby : {groupby}
    type_order : {type_order}
    key : {plot_key}
    **plot_options :
        Pass to :func:`milkviz.dot` and :func:`milkviz.anno_clustermap`.

    """
    pdata = get_result(data, key)
    pdata = pd.pivot_table(
        pdata, columns="cell_type", index=pdata.index.names[1::], values="pattern"
    )
    if type_order is None:
        type_order = natsorted(pdata.columns)
    pdata = pdata[type_order]

    if use == "dot":
        pattern = {}
        for t, arr in pdata.iteritems():
            pattern[t] = {0: 0, 1: 0, 2: 0, 3: 0, **Counter(arr)}
        pdata = pd.DataFrame(pattern).T[[0, 1, 2, 3]]
        colors = np.repeat(
            [["#FFC408", "#c54a52", "#4a89b9", "#5a539d"]], len(pdata), axis=0
        )
        return dot(
            dot_size=pdata.to_numpy(dtype=int),
            dot_hue=colors,
            xticklabels=["No Cell", "Random", "Regular", "Cluster"],
            yticklabels=pdata.index,
            **plot_options,
        )
    else:
        pdata = pdata.rename_axis(columns={"cell_type": "Cell Type"})
        plot_kw = dict(
            categorical_cbar=["No Cell", "Random", "Regular", "Cluster"],
            heat_cmap=ListedColormap(["#FFC408", "#c54a52", "#4a89b9", "#5a539d"]),
            col_legend_split=False,
            cbar_title="Pattern",
            vmin=0,
            vmax=3,
            col_cluster=False,
        )
        plot_kw = {**plot_kw, **plot_options}
        return anno_clustermap(
            pdata,
            col_colors="Cell Type",
            row_colors=groupby,
            **plot_kw
        )
