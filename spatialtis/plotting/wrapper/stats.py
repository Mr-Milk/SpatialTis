from collections import OrderedDict
from itertools import combinations_with_replacement
from typing import Optional, Sequence

import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import heatmap, stacked_bar, tri_dotplot, violin_plot
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def cell_components(
    adata: AnnData,
    groupby: Sequence[str],
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) Plotting function for cell components

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        **kwargs: Pass to `Stacked-bar plot <plotting.html#spatialtis.plotting.stacked_bar>`_

    """
    if key is None:
        key = CONFIG.cell_components_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[selected_types]

    p = stacked_bar(df, groupby, **kwargs)

    return p


@reuse_docstring()
def cell_density(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) plotting function for cell density

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        **kwargs: Pass to `Violin plot <plotting.html#spatialtis.plotting.violin_plot>`_

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

    df = df.iloc[df.index.sortlevel(groupby)[1], :]

    levels_order = groupby.copy()
    for name in df.index.names:
        if name not in levels_order:
            levels_order.append(name)

    df = df.reorder_levels(levels_order)

    p = violin_plot(df, groupby, "density", **kwargs)

    return p


@reuse_docstring()
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
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        use: Options are "dot" and "heatmap"
        key: {key}
        **kwargs: Pass to `Tri-Dotplot <plotting.html#spatialtis.plotting.tri_dotplot>`_ or
                  `Heatmap <plotting.html#spatialtis.plotting.heatmap>`_

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
        plot_kwargs = dict(legend_title="ROI")

        for k, v in kwargs.items():
            plot_kwargs[k] = v
        p = tri_dotplot(counts, labels=labels, **plot_kwargs)
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


@reuse_docstring()
def cell_morphology(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """(bokeh) plotting function for cell morphology

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        **kwargs: Pass to `Violin plot <plotting.html#spatialtis.plotting.violin_plot>`_

    """
    if key is None:
        key = CONFIG.cell_morphology_key

    df = adata_uns2df(adata, key)

    if selected_types is not None:
        df = df[df.index.get_level_values("type").isin(selected_types)]

    if groupby is not None:
        if "type" not in groupby:
            groupby = ["type"] + list(groupby)
        else:
            groupby = list(groupby)
    else:
        groupby = ["type"]

    df = df.iloc[df.index.sortlevel(groupby)[1], :]

    levels_order = groupby.copy()
    for name in df.index.names:
        if name not in levels_order:
            levels_order.append(name)

    df = df.reorder_levels(levels_order)

    p = violin_plot(df, groupby, "value", **kwargs)

    return p
