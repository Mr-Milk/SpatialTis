from itertools import combinations
from typing import Optional, Sequence, Mapping

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import chisquare

from spatialtis.config import CONFIG
from spatialtis.utils import adata_uns2df

from ._bar_plot import stacked_bar
from ._heatmap_sns import heatmap
from ._stacked_kde_sns import stacked_kde
from ._violin_plot import violin_plot
from ._cell_cell_interaction import cc_interactions
from ._grouped_pie import grouped_pie
from ._community_graph import graph_plot


def cell_components(
        adata: AnnData,
        groupby: Sequence[str],
        selected_types: Optional[Sequence] = None,
        key: str = "cell_components",
        **kwargs
):
    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    p = stacked_bar(df, groupby, **kwargs)

    return p


def cell_density(
        adata: AnnData,
        groupby: Sequence[str],
        selected_types: Optional[Sequence] = None,
        key: str = "cell_density",
        **kwargs
):
    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    df = pd.DataFrame(df.stack(), columns=["density"])
    groupby = list(groupby) + ["type"]
    p = violin_plot(df, groupby, "density", **kwargs)

    return p


def cell_co_occurrence(
        adata: AnnData,
        groupby: Sequence[str],
        selected_types: Optional[Sequence] = None,
        key: str = "cell_co_occurrence",
        pval: float = 0.01,
        **kwargs
):
    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    tdf = df.astype(int).groupby(level=groupby)
    X = []
    for n, g in tdf:
        c = (g.sum() / len(g)).fillna(0)
        data = []
        for comb in combinations(c, 2):
            if (comb[0] == 0) | (comb[1] == 0):
                p = 0
            else:
                p = chisquare(comb).pvalue
            data.append(p)
        X.append(data)
    pdf = (pd.DataFrame(X) > (1 - pval)).astype(int)
    pdf.index = pd.MultiIndex.from_frame(
        df.index.to_frame(index=False)[groupby].drop_duplicates()
    )
    pdf.columns = pd.MultiIndex.from_arrays(
        np.asarray([i for i in combinations(df.columns, 2)]).T,
        names=["Cell type1", "Cell type2"],
    )

    plot_kwargs = dict(
        row_colors=groupby,
        col_colors=["Cell type1", "Cell type2"],
        colorbar_type="categorical",
        categorical_colorbar_text=["Absent", "Presence"],
        col_colors_legend_bbox=(1.05, 0.5),
        row_colors_legend_bbox=(-0.25, 0.5),
        colorbar_bbox=(-0.25, 0.15),
        row_cluster=None,
        col_cluster=True,
    )
    # allow user to overwrite the default plot config
    for k, v in kwargs.items():
        plot_kwargs[k] = v

    p = heatmap(pdf, **plot_kwargs)
    return p


def cell_morphology(
        adata: AnnData,
        row: Optional[str] = None,
        col: Optional[str] = None,
        selected_types: Optional[str] = None,
        *,
        key: str = "cell_morphology",
        **kwargs
):
    df = adata_uns2df(adata, key)
    type_col = CONFIG.CELL_TYPE_COL

    p = stacked_kde(df.xs(selected_types, level=type_col), row=row, col=col, **kwargs)
    return p


def neighborhood_analysis(adata: AnnData,
                          groupby: Sequence[str],
                          key: str = "neighborhood_analysis",
                          method: str = "graph",  # graph, heatmap
                          **kwargs
                          ):
    df = adata_uns2df(adata, key)

    if method == 'graph':
        combs = df.columns
        new_cols = list()
        for comb in combs:
            add_numbers = ((comb[0] + '-1'), (comb[1] + '-2'))
            new_cols.append(add_numbers)

        df.columns = new_cols

        p = cc_interactions(df, {-1: "Avoidance", 1: "Association"}, order=[-1, 1], **kwargs)
    elif method == 'heatmap':
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["Cell type1", "Cell type2"],
            palette=["#2f71ab", "#f7f7f7", "#ba262b"],
            colorbar_type="categorical",
            categorical_colorbar_text=["Avoidance", "Association"],
            col_colors_legend_bbox=(1.05, 0.5),
            row_colors_legend_bbox=(-0.25, 0.5),
            colorbar_bbox=(-0.25, 0.15),
            row_cluster=None,
            col_cluster=True,
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = heatmap(df, **plot_kwargs)
    return p


def spatial_enrichment_analysis(
        adata: AnnData,
        groupby: Sequence[str],
        key: str = "spatial_enrichment_analysis",
        **kwargs
):
    df = adata_uns2df(adata, key)

    plot_kwargs = dict(
        row_colors=groupby,
        col_colors=["Cell type1", "Cell type2"],
        col_colors_legend_bbox=(1.05, 0.5),
        row_colors_legend_bbox=(-0.25, 0.15),
        row_cluster=None,
        col_cluster=True,
    )
    # allow user to overwrite the default plot config
    for k, v in kwargs.items():
        plot_kwargs[k] = v

    p = heatmap(df, **plot_kwargs)
    return p


def spatial_distribution(adata: AnnData,
                         groupby: Sequence[str],
                         key: str = "spatial_distribution",
                         method: str = 'pie',  # pie, heatmap
                         **kwargs
                         ):
    df = adata_uns2df(adata, key)

    if method == 'pie':
        order = [1, 2, 3, 0]
        mapper = {0: "No cell", 1: "Random", 2: "Regular", 3: "Cluster"}
        p = grouped_pie(df, mapper, order=order, **kwargs)
    elif method == 'heatmap':
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["Cell type"],
            palette=["#fffec6", "#c54a52", "#4a89b9", "#5a539d"],
            colorbar_type="categorical",
            categorical_colorbar_text=["Blank", "Random", "Regular", "Clumped"],
            col_colors_legend_bbox=(1.05, 0.5),
            row_colors_legend_bbox=(-0.25, 0.5),
            row_cluster=None,
            col_cluster=True,
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = heatmap(df, **plot_kwargs)
    return p


def spatial_heterogeneity(
        adata: AnnData,
        groupby: Sequence[str],
        key: str = "spatial_heterogeneity",
        metric: str = "heterogeneity",
        **kwargs
):
    df = adata_uns2df(adata, key)

    p = violin_plot(df, groupby, metric, **kwargs)

    return p


def cell_type_graph(
        adata: AnnData,
        query: Mapping,
        type_col: Optional[str] = None,
        neighbors_col: Optional[str] = None,
        centroid_col: Optional[str] = None,
        **kwargs,
):

    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if neighbors_col is None:
        neighbors_col = CONFIG.NEIGHBORS_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    p = graph_plot(df,
                   node_col=centroid_col,
                   node_category_col=type_col,
                   edge_col=neighbors_col,
                   **kwargs)
    return p


def cell_communities_graph(
        adata: AnnData,
        query: Mapping,
        type_col: Optional[str] = None,
        community_col: Optional[str] = None,
        neighbors_col: Optional[str] = None,
        centroid_col: Optional[str] = None,
        **kwargs,
):
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if community_col is None:
        community_col = CONFIG.COMMUNITY_COL
    if neighbors_col is None:
        neighbors_col = CONFIG.NEIGHBORS_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    p = graph_plot(df,
                   node_col=centroid_col,
                   node_info_col=type_col,
                   edge_col=neighbors_col,
                   edge_category_col=community_col,
                   **kwargs)
    return p
