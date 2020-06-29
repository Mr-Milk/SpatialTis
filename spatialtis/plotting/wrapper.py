from collections import Counter, OrderedDict
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
    use: str = "dot",  # dot, heatmap
    key: Optional[str] = None,
    **kwargs,
):
    """(matplotlib) plotting function for cell co-occurrence

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        use: "dot" or "heatmap"
        key: which key to read the data
        **kwargs: pass to plotting.tri_dotplot (use="dot") or plotting.heatmap (use="heatmap")

    Returns:

    """
    if key is None:
        key = CONFIG.cell_co_occurrence_key

    df = adata_uns2df(adata, key)

    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if selected_types is not None:
        df = df[selected_types]
        df = df.iloc[:, df.columns.get_level_values("Cell type2").isin(selected_types)]

    if use == "dot":
        ndf = df.reset_index(drop=True).sum().reset_index()
        labels = pd.unique(ndf["Cell type1"])
        if selected_types is not None:
            labels = selected_types
        combs = [i for i in combinations_with_replacement(labels, 2)]
        counts = OrderedDict((k, []) for k in labels)

        for comb in combs:
            v = ndf[(ndf["Cell type1"] == comb[0]) & (ndf["Cell type2"] == comb[1])]
            if len(v) == 0:
                v = ndf[(ndf["Cell type1"] == comb[1]) & (ndf["Cell type2"] == comb[0])]
            v = v[0].tolist()[0]
            counts[comb[0]].append(v)

        counts = list(counts.values())
        p = tri_dotplot(counts, labels, **kwargs)
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
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "dot_matrix",  # graph, heatmap, dot_matrix
    **kwargs,
):
    """("dot_matrix", "heatmap": matplotlib, "graph": pyechart) plotting function for neighborhood analysis

    Args:
        adata: anndata object
        groupby: how to group your data in plot, only work for use="heatmap"
        selected_types:
        key: which key to read the data
        use: "dot_matrix", "heatmap", "graph"
        **kwargs: pass to plotting.dot_matrix; plotting.heatmap; plotting.cc_interaction

    """
    if key is None:
        key = CONFIG.neighborhood_analysis_key

    df, params = adata_uns2df(adata, key, params=True)
    order = params["order"]

    if selected_types is not None:
        combs = [i for i in combinations_with_replacement(selected_types, 2)]
        cols = df.columns.tolist()
        colnames = df.columns.names
        for i, comb in enumerate(cols):
            if comb not in combs:
                cols[i] = comb[::-1]
        df.columns = pd.MultiIndex.from_tuples(cols)
        df.columns.names = colnames
        df = df[combs]

    if use == "graph":
        p = cc_interactions(df, {-1: "Avoidance", 1: "Association"}, **kwargs)
    elif use == "dot_matrix":
        if not order:
            raise ValueError(
                "Unordered cell-interaction couldn't be visualized using dot matrix plot"
            )
        try:
            cc = adata_uns2df(adata, CONFIG.cell_components_key)
        except KeyError:
            raise Exception(
                "Please run cell_components before plotting the dot matrix plot"
            )
        counts = (
            list()
        )  # [cell 1, cell 2, interaction type, No. of sig (both type), % of sign type, p)
        for label, data in df.items():
            c_data = sorted(
                [(k, v) for k, v in {**{0: 0, 1: 0, -1: 0}, **Counter(data)}.items()],
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
            if all_sign == 0:
                counts.append([label[0], label[1], 0, 0, p])
            else:
                if c_data[0][0] == 0:
                    if c_data[1][0] == 1:
                        counts.append(
                            [label[0], label[1], all_sign, c_data[1][1] / all_sign, p]
                        )
                    else:
                        counts.append(
                            [
                                label[0],
                                label[1],
                                all_sign,
                                1 - c_data[1][1] / all_sign,
                                p,
                            ]
                        )
                elif c_data[0][0] == 1:
                    counts.append(
                        [label[0], label[1], all_sign, c_data[0][1] / all_sign, p]
                    )
                else:
                    counts.append(
                        [label[0], label[1], all_sign, 1 - c_data[0][1] / all_sign, p]
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

    elif use == "heatmap":
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
        unique_values = np.unique(df.to_numpy())
        if 1 not in unique_values:
            plot_kwargs["palette"] = ["#2f71ab", "#f7f7f7"]
            plot_kwargs["categorical_colorbar_text"] = ["Avoidance", "no-sign"]
        if -1 not in unique_values:
            plot_kwargs["palette"] = ["#f7f7f7", "#ba262b"]
            plot_kwargs["categorical_colorbar_text"] = ["no-sign", "Association"]
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
    use: str = "dot",  # heatmap
    **kwargs,
):
    """(matplotlib) plotting function for spatial distribution

    Args:
        adata: anndata object
        groupby: how to group your data in plot
        selected_types: select interested types
        key: which key to read the data
        use: "dot" or "heatmap"
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
    if use == "dot":
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
        colors = np.array(["#FFC408", "#c54a52", "#4a89b9", "#5a539d"] * len(tb))

        plot_kwargs = dict(annotated=False, color=colors, alpha=1,)

        for k, v in kwargs:
            plot_kwargs[k] = v

        p = dotplot(tb, y="Cell", x="Pattern", **plot_kwargs)
        # p = dotplot(tb.T, x="Cell", y="Pattern", annotated=False)
    elif use == "heatmap":
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
        raise ValueError("Support plotting methods are 'dot' and 'heatmap'.")
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


def cell_neighbors(
    adata: AnnData,
    query: Mapping,
    use: str = "interactive",
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
    df = df.reset_index()

    nodes = [eval(n) for n in df[centroid_key]]
    nodes_types = [str(n) for n in df[type_key]]

    neighs = [eval(n) for n in df[neighbors_key]]
    edges = []
    for i, n in zip(df.index, neighs):
        for x in n:
            edges.append((i, x))

    if use == "interactive":
        p = graph_plot_interactive(nodes, edges, nodes_types=nodes_types, node_size=3)
    else:
        p = graph_plot(nodes, edges, nodes_types=nodes_types, node_size=4)
    return p


def cell_communities(
    adata: AnnData,
    query: Mapping,
    min_cell: int = 10,
    use: str = "interactive",
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
        min_cell: only show communities that have certain number of cells
        use: "interactive" or "static", for big ROI, "interactive" is much faster using WebGL
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

    nodes_types = df[community_key].tolist()
    commus = []
    for commu, count in Counter(nodes_types).items():
        if count >= min_cell:
            commus.append(commu)

    df = df.reset_index(drop=True)
    xdf = df[df[community_key].isin(commus)]
    xdf = xdf.reset_index()

    nodes = [eval(n) for n in xdf[centroid_key]]
    neighs = [eval(n) for n in xdf[neighbors_key]]
    nodes_types = xdf[community_key]

    edges = []
    edges_types = []
    for i, n in zip(xdf.index, neighs):
        for x in n:
            new_x = xdf[xdf["index"] == x].index
            if len(new_x) == 1:
                new_x = new_x[0]
                if nodes_types[i] == nodes_types[new_x]:
                    edges.append((i, new_x))
                    edges_types.append(nodes_types[i])

    """
    nodes = []
    edges = []
    edges_types = []
    for i, d in xdf.iterrows():
        nodes.append(eval(d[centroid_key]))
        neighs = eval(d[neighbors_key])
        t = d[community_key]

        for n in neighs:
            info = xdf[xdf['index'] == n][community_key]
            if len(info) == 1:
                if t == info.tolist()[0]:
                    edges.append((i, info.index[0]))
                    edges_types.append(t)
"""
    if use == "interactive":
        p = graph_plot_interactive(nodes, edges, edges_types=edges_types, **kwargs)
    else:
        p = graph_plot(nodes, edges, edges_types=edges_types, **kwargs)

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
