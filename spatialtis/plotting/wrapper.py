from collections import Counter
from itertools import combinations_with_replacement
from typing import Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import chisquare, pearsonr

from spatialtis.config import CONFIG
from spatialtis.plotting.palette import get_colors
from spatialtis.utils import adata_uns2df

from ._bar_plot import stacked_bar
from ._cell_cell_interaction import cc_interactions
from ._community_graph import graph_plot, graph_plot_interactive
from ._dot_matrixplot import dot_matrix
from ._dotplot import dotplot
from ._heatmap_sns import heatmap
from ._sankey import sankey
from ._stacked_kde_sns import stacked_kde
from ._triangle_dotplot import tri_dotplot
from ._violin_plot import violin_plot


def cell_components(
    adata: AnnData,
    groupby: Sequence[str],
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) plotting function for cell components

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        key: which key to read the data
        **kwargs: pass to plotting.stacked_bar

    """
    if key is None:
        key = CONFIG.cell_components_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    p = stacked_bar(df, groupby, **kwargs)

    return p


def cell_density(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) plotting function for cell density

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        key: which key to read the data
        **kwargs: pass to plotting.violin_plot

    """
    if key is None:
        key = CONFIG.cell_density_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    df = pd.DataFrame(df.stack(), columns=["density"])
    if groupby is not None:
        if "type" not in groupby:
            groupby = ["type"] + list(groupby)
        else:
            groupby = list(groupby)
    else:
        groupby = ["type"]
    p = violin_plot(df, groupby, "density", **kwargs)

    return p


def cell_co_occurrence(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    method: str = "dot",  # dot, heatmap
    key: Optional[str] = None,
    **kwargs,
):
    """(matplotlib) plotting function for cell co-occurrence

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        method: "dot" or "heatmap"
        key: which key to read the data
        **kwargs: pass to plotting.tro_dotplot (method="dot") or plotting.heatmap (method="heatmap")

    Returns:

    """
    if key is None:
        key = CONFIG.cell_co_occurrence_key

    df = adata_uns2df(adata, key)

    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if selected_types is not None:
        df = df[selected_types]

    if method == "dot":
        ndf = df.reset_index(drop=True).sum().reset_index()
        labels = pd.unique(ndf["Cell type1"])
        sorted_labels = []
        if selected_types is not None:
            for t in selected_types:
                if t in labels:
                    sorted_labels.append(t)
        else:
            sorted_labels = labels

        counts = list()
        for n, g in ndf.groupby("Cell type1"):
            arr = list(g[0].to_numpy())
            counts.append(arr)

        counts = sorted(counts, key=len, reverse=True)
        p = tri_dotplot(counts, sorted_labels, **kwargs)
    else:
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

        p = heatmap(df, **plot_kwargs)

    return p


def cell_morphology(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) plotting function for cell morphology

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        key: which key to read the data
        **kwargs: pass to plotting.violin_plot

    """
    if key is None:
        key = CONFIG.cell_morphology_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[df["type"].isin(selected_types)]

    if groupby is not None:
        df = df.set_index(groupby + ["type", "id"], drop=True)
        if "type" not in groupby:
            groupby = ["type"] + list(groupby)
        else:
            groupby = list(groupby)
    else:
        df = df.set_index(["type", "id"], drop=True)
        groupby = ["type"]
    p = violin_plot(df, groupby, "value", **kwargs)

    return p


def neighborhood_analysis(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    key: Optional[str] = None,
    method: str = "dot_matrix",  # graph, heatmap, dot_matrix
    **kwargs,
):
    """("dot_matrix", "heatmap": matplotlib, "graph": pyechart) plotting function for neighborhood analysis

    Args:
        adata: anndata object
        groupby: how to group your data in plot, only work for method="heatmap"
        key: which key to read the data
        method: "dot_matrix", "heatmap", "graph"
        **kwargs: pass to plotting.dot_matrix; plotting.heatmap; plotting.cc_interaction

    """
    if key is None:
        key = CONFIG.neighborhood_analysis_key

    df = adata_uns2df(adata, key)

    if method == "graph":
        p = cc_interactions(df, {-1: "Avoidance", 1: "Association"}, **kwargs)
    elif method == "dot_matrix":
        try:
            cc = adata_uns2df(adata, CONFIG.cell_components_key)
        except KeyError:
            raise Exception(
                "Please run cell_components before plotting the dot matrix plot"
            )

        counts = (
            list()
        )  # [cell 1, cell 2, interaction type, No. of sig (both type), % of sign type, p)
        for i, (label, data) in enumerate(df.items()):
            c_data = sorted(
                [(k, v) for k, v in Counter(data).items()],
                key=lambda x: x[1],
                reverse=True,
            )

            # count the No. of sign for both type
            all_sign = 0
            for (t, n) in c_data:
                if t != 0:
                    all_sign += n

            # count pearson correlation between cell components
            p = pearsonr(cc[label[0]], cc[label[1]])[0]

            if c_data[0][0] == 0:
                counts.append([label[0], label[1], all_sign, 0, p])
            elif c_data[0][0] == 1:
                counts.append(
                    [label[0], label[1], all_sign, 1 - c_data[0][1] / all_sign, p]
                )
            else:
                counts.append(
                    [label[0], label[1], all_sign, c_data[0][1] / all_sign, p]
                )
        counts = pd.DataFrame(counts, columns=["type1", "type2", "all", "%", "p"])
        matrix = counts.pivot(index="type1", columns="type2", values="p")
        xlabels = matrix.columns
        ylabels = matrix.index
        matrix = matrix.to_numpy()
        dot_color = counts.pivot(index="type1", columns="type2", values="%").to_numpy()
        dot_size = counts.pivot(index="type1", columns="type2", values="all").to_numpy()

        plot_kwargs = dict(
            size_legend_title="Sign 'ROI Counts",
            color_legend_title="% of interaction",
            cbar_legend_title="Pearson\nCorrelation",
            cbar_mapper={-1: "-1", 1: "1"},
            color_legend_text=["Association", "", "Avoidance"],
        )
        for k, v in kwargs.items():
            plot_kwargs[k] = v
        p = dot_matrix(matrix, dot_color, dot_size, xlabels, ylabels, **plot_kwargs)

    elif method == "heatmap":
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["Cell type1", "Cell type2"],
            palette=["#2f71ab", "#f7f7f7", "#ba262b"],
            colorbar_type="categorical",
            categorical_colorbar_text=["Avoidance", "Association"],
            col_colors_legend_bbox=(1.05, 0.5),
            row_colors_legend_bbox=(-0.25, 0.5),
            colorbar_bbox=(-0.25, 0.15),
            row_cluster=True,
            col_cluster=None,
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = heatmap(df, **plot_kwargs)
    return p


def spatial_enrichment_analysis(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(matplotlib) plotting function for spatial enrichment analysis

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        key: which key to read the data
        **kwargs: pass to plotting.heatmap

    """
    if key is None:
        key = CONFIG.spatial_enrichment_analysis_key

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


def spatial_distribution(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    method: str = "dot",  # heatmap
    **kwargs,
):
    """(matplotlib) plotting function for spatial distribution

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        key: which key to read the data
        method: "dot" or "heatmap"
        **kwargs: pass to plotting.dotplot or plotting.heatmap

    """
    if key is None:
        key = CONFIG.spatial_distribution_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    """ Overlay text label, maybe support in future
    if method == "pie":
        order = [1, 2, 3, 0]
        mapper = {0: "No cell", 1: "Random", 2: "Regular", 3: "Cluster"}
        p = grouped_pie(df, mapper, order=order, **kwargs)
    """
    if method == "dot":
        names = []
        counts = []
        for n, col in df.iteritems():
            names.append(n)
            counts.append(Counter(col))
        counts.append({1: 0, 2: 0, 3: 0, 0: 0})
        tb = pd.DataFrame(counts)
        tb.drop(tb.tail(1).index, inplace=True)
        tb.fillna(0, inplace=True)
        tb.rename(
            columns={0: "No Cell", 1: "Random", 2: "Regular", 3: "Cluster"},
            inplace=True,
        )
        tb.index = names
        tb.index.name = "Cell"
        tb.columns.name = "Pattern"
        tb = tb[["No Cell", "Random", "Regular", "Cluster"]]
        p = dotplot(tb, y="Cell", x="Pattern", annotated=False, **kwargs)
        # p = dotplot(tb.T, x="Cell", y="Pattern", annotated=False)
    elif method == "heatmap":
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["Cell type"],
            palette=["#fffec6", "#c54a52", "#4a89b9", "#5a539d"],
            colorbar_type="categorical",
            categorical_colorbar_text=["No Cell", "Random", "Regular", "Cluster"],
            col_colors_legend_bbox=(1.05, 0.5),
            row_colors_legend_bbox=(-0.25, 0.5),
            row_cluster=None,
            col_cluster=True,
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = heatmap(df, **plot_kwargs)
    else:
        raise ValueError("Support methods are 'dot' and 'heatmap'.")
    return p


def spatial_heterogeneity(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    key: Optional[str] = None,
    metric: str = "heterogeneity",
    **kwargs,
):
    """(bokeh) plotting function for cell morphology

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        key: which key to read the data
        metric: "heterogeneity" or "KL", "KL" only available if you use shannon entropy
        **kwargs: pass to plotting.violin_plot or plotting.stacked_bar

    """
    if key is None:
        key = CONFIG.spatial_heterogeneity_key

    if metric not in ["heterogeneity", "KL"]:
        raise ValueError("Available options for metric are 'heterogeneity' and 'KL'.")

    df = adata_uns2df(adata, key)
    inames = df.index.names

    if len(inames) == 2:
        p = stacked_bar(df, [inames[0]], percentage=False, sort_type=metric)
    else:
        p = violin_plot(df, groupby, metric, **kwargs)

    return p


def cell_type_graph(
    adata: AnnData,
    query: Mapping,
    type_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    **kwargs,
):
    """(pyecharts) visualize cell type in ROI

    Args:
        adata:
        query:
        type_key: key to cell type
        neighbors_key: key to neighbors info
        centroid_key: key to cell centroid
        **kwargs: pass to plotting.graph_plot_interactive

    Returns:

    """
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.neighbors_key
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    p = graph_plot_interactive(
        df,
        node_col=centroid_key,
        node_category_col=type_key,
        edge_col=neighbors_key,
        **kwargs,
    )
    return p


def cell_communities_graph(
    adata: AnnData,
    query: Mapping,
    method: str = "interactive",
    type_key: Optional[str] = None,
    community_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    **kwargs,
):
    """("interactive": pyecharts, "static": matplotlib) Visualize cell spatial communiteis

    Args:
        adata: anndata
        query: a dict use to select which ROI to display
        method: "interactive" or "static", for big ROI, "interactive" is much faster using WebGL
        type_key: key to cell type
        community_key: key to community
        neighbors_key: key to neighbors
        centroid_key: key to cell centroid
        **kwargs: pass to plotting.graph_plot_interactive or plotting.graph_plot


    """
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if community_key is None:
        community_key = CONFIG.community_key
    if neighbors_key is None:
        neighbors_key = CONFIG.neighbors_key
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    if neighbors_key not in adata.obs.keys():
        raise KeyError(
            "Neighbor key not found, export neighbors or specific your own key."
        )
    if community_key not in adata.obs.keys():
        raise KeyError(
            "Community key not found, run community or specific your own key."
        )

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))

    if method == "interactive":
        p = graph_plot_interactive(
            df,
            node_col=centroid_key,
            node_info_col=type_key,
            edge_col=neighbors_key,
            edge_category_col=community_key,
            **kwargs,
        )

    else:
        nodes = [eval(n) for n in df[centroid_key]]
        neighs = df[neighbors_key]
        nodes_types = list(df[community_key])

        edges = []
        edges_types = []
        for i, n in enumerate(neighs):
            for x in n:
                if nodes_types[i] == nodes_types[x]:
                    edges_types.append(nodes_types[i])
                    edges.append((i, x))

        p = graph_plot(nodes, edges, edges_types=edges_types)

    return p


def exp_neighcells(
    adata: AnnData,
    key: Optional[str] = None,
    palette: Optional[Sequence] = None,
    **kwargs,
):
    """(pyecharts) plotting function for expression influenced by neighbor cells

    Args:
        adata: anndata object
        key: key to read data
        palette: config the color
        **kwargs: pass to plotting.sankey

    """
    if key is None:
        key = CONFIG.exp_neighcell_key

    df = adata_uns2df(adata, key)

    cell_types = pd.unique(df.iloc[:, [0, 1]].values.flatten())
    gene_types = pd.unique(df.iloc[:, [2]].values.flatten())

    if palette is None:
        palette = ["Set3", "Spectral"]

    colors = get_colors(len(cell_types) + len(gene_types), palette)
    cell_colormap = dict(zip(cell_types, colors[0 : len(cell_types)]))
    gene_colormap = dict(
        zip(gene_types, colors[len(cell_types) : (len(cell_types) + len(gene_types))])
    )

    nodes = []
    nodes_colors = []
    links = []
    for c in pd.unique(df.iloc[:, [0]].values.flatten()):
        nodes.append(c)
        nodes_colors.append(cell_colormap[c])
    for c in pd.unique(df.iloc[:, [1]].values.flatten()):
        nodes.append(c + " ")
        nodes_colors.append(cell_colormap[c])
    for g in gene_types:
        nodes.append(g)
        nodes_colors.append(gene_colormap[g])
    for i, row in df.iterrows():
        links.append((row[0], row[1] + " ", 0.1))
        links.append((row[1] + " ", row[2], row[3]))

    p = sankey(nodes, nodes_colors, links, **kwargs)

    # return nodes, nodes_colors, links
    return p
