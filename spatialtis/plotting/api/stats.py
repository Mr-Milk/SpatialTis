from typing import List, Optional, Sequence

import numpy as np
from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import (
    dot_plot,
    heatmap,
    stacked_bar_interactive,
    stacked_bar_static,
    violin_static,
)
from spatialtis.typing import Array
from spatialtis.utils import doc, get_result


@doc
def cell_components(
    data: AnnData,
    groupby: List[str],
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "static",
    **kwargs,
):
    """Visualization for cell components

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        use: "static" or "interactive" (Default: "static")
        **kwargs: Pass to :class:`spatialtis.plotting.base.stacked_bar_static` or
            :class:`spatialtis.plotting.base.stacked_bar_interactive`

    """
    if key is None:
        key = ANALYSIS["cell_components"].last_used_key

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]

    if selected_types is not None:
        df = df[exp_obs + selected_types]
    stacked_types = list(df.columns)
    for i in exp_obs:
        stacked_types.remove(i)

    if use == "interactive":
        return stacked_bar_interactive(
            df, groupby, stacked_types, **kwargs, saved_name="cell_components"
        )
    else:
        return stacked_bar_static(
            df, groupby, stacked_types, **kwargs, saved_name="cell_components"
        )


@doc
def cell_density(
    data: AnnData,
    groupby: List[str],
    selected_types: Optional[List] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """Visualization for cell density

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        **kwargs: Pass to :class:`spatialtis.plotting.base.violin_static`

    """
    gl = len(groupby)
    if gl > 2:
        raise ValueError("Only support 2 levels depth categorical data for matplotlib")
    if key is None:
        key = ANALYSIS["cell_density"].last_used_key

    df = get_result(data, key)

    if selected_types is not None:
        df = df[df["type"].isin(selected_types)]

    return violin_static(
        df, groupby, target="value", hue="type", **kwargs, saved_name="cell_density"
    )


@doc
def cell_co_occurrence(
    data: AnnData,
    groupby: Optional[Array] = None,
    selected_types: Optional[Sequence] = None,
    use: str = "dot",  # dot, heatmap
    key: Optional[str] = None,
    **kwargs,
):
    """Visualization for cell co-occurrence

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        use: "dot" or "heatmap" (Default: "dot")
        key: {key}
        **kwargs: Pass to :class:`spatialtis.plotting.base.dot_plot` or
            :class:`spatialtis.plotting.base.heatmap`

    """
    if key is None:
        key = ANALYSIS["cell_co_occurrence"].last_used_key

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]

    if groupby is None:
        groupby = exp_obs
    if selected_types is not None:
        df = df[df["type1"].isin(selected_types) & df["type2"].isin(selected_types)]

    if use == "dot":
        ndf = df.pivot_table(
            index="type1", columns="type2", values="co_occur", aggfunc=np.sum
        )
        order = ndf.count().sort_values().index.tolist()

        plot_kwargs = dict(legend_title="ROI")
        for k, v in kwargs.items():
            plot_kwargs[k] = v
        p = dot_plot(
            ndf.loc[order, order[::-1]], **plot_kwargs, saved_name="cell_co_occurrence"
        )
    else:
        ndf = df.pivot_table(
            index=exp_obs, columns=["type1", "type2"], values="co_occur", aggfunc=np.sum
        )
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["type1", "type2"],
            categorical_colorbar=["Absent", "Presence"],
            clustermap_kwargs=dict(row_cluster=None, col_cluster=True,),
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = heatmap(ndf, **plot_kwargs, saved_name="cell_co_occurrence")

    return p


@doc
def cell_morphology(
    data: AnnData,
    groupby: List[str],
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    **kwargs,
):
    """Visualization for cell morphology

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        **kwargs: Pass to :class:`spatialtis.plotting.base.violin_static`

    """
    gl = len(groupby)
    if gl > 2:
        raise ValueError("Only support 2 levels depth categorical data for matplotlib")
    if key is None:
        key = ANALYSIS["cell_morphology"].last_used_key

    df = get_result(data, key)

    if selected_types is not None:
        df = df[df["type"].isin(selected_types)]

    return violin_static(
        df, groupby, target="value", hue="type", **kwargs, saved_name="cell_morphology",
    )
