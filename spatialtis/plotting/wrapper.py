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
from ._grouped_pie import grouped_pie
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
    if key is None:
        key = CONFIG.cell_density_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    df = pd.DataFrame(df.stack(), columns=["density"])
    if groupby is not None:
        groupby = ["type"] + list(groupby)
    else:
        groupby = ["type"]
    p = violin_plot(df, groupby, "density", **kwargs)

    return p


def cell_co_occurrence(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    compare: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    method: str = "dot",  # dot, heatmap
    key: Optional[str] = None,
    pval: float = 0.01,
    **kwargs,
):
    if key is None:
        key = CONFIG.cell_co_occurrence_key

    df = adata_uns2df(adata, key)

    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if compare is None:
        compare = CONFIG.EXP_OBS[0]
    if selected_types is not None:
        df = df[selected_types]

    tdf = df.astype(int).groupby(level=compare)
    X = []
    for n, g in tdf:
        c = (g.sum() / len(g)).fillna(0)
        data = []
        for comb in combinations_with_replacement(c, 2):
            if (comb[0] == 0) | (comb[1] == 0):
                p = 0
            else:
                p = chisquare(comb).pvalue
            data.append(p)
        X.append(data)
    pdf = (pd.DataFrame(X) > (1 - pval)).astype(int)
    pdf.index = df.index.to_frame(index=False)[compare].drop_duplicates()
    pdf.columns = pd.MultiIndex.from_arrays(
        np.asarray([i for i in combinations_with_replacement(df.columns, 2)]).T,
        names=["Cell type1", "Cell type2"],
    )

    if method == "dot":
        ndf = pdf.reset_index(drop=True).sum().reset_index()
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

        p = heatmap(pdf, **plot_kwargs)

    return p


def cell_morphology(
    adata: AnnData,
    row: Optional[str] = None,
    col: Optional[str] = None,
    selected_types: Optional[str] = None,
    *,
    key: Optional[str] = None,
    **kwargs,
):
    if key is None:
        key = CONFIG.cell_morphology_key

    df = adata_uns2df(adata, key)
    type_col = CONFIG.CELL_TYPE_KEY

    p = stacked_kde(df.xs(selected_types, level=type_col), row=row, col=col, **kwargs)
    return p


def neighborhood_analysis(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    key: Optional[str] = None,
    method: str = "graph",  # graph, heatmap, dot_matrix
    order: bool = False,
    **kwargs,
):
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
        p = dot_matrix(
            matrix,
            dot_color,
            dot_size,
            xlabels,
            ylabels,
            size_legend_title="Sign 'ROI Counts",
            color_legend_title="% of interaction",
            cbar_legend_title="Pearson\nCorrelation",
            cbar_mapper={-1: "-1", 1: "1"},
            color_legend_text=["Association", "", "Avoidance"],
        )

    elif method == "heatmap":

        return df

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
    method: str = "dot",  # pie, heatmap
    **kwargs,
):
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
        p = dotplot(tb, y="Cell", x="Pattern", annotated=False)
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
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.NEIGHBORS_KEY
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
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if community_key is None:
        community_key = CONFIG.COMMUNITY_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.NEIGHBORS_KEY
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
