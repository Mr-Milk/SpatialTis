from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import (
    dot_matrix,
    graph_layout_interactive,
    graph_layout_static,
    heatmap,
)
from spatialtis.utils import doc, get_result


@doc
def spatial_enrichment_analysis(
    data: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_markers: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "dot",  # heatmap, graph
    **kwargs,
):
    """Visualization for spatial enrichment analysis

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_markers: {selected_markers}
        key: {key}
        use: "dot", "heatmap", "graph_static" or "graph_interactive" (Default: "dot")
        **kwargs: Pass to :class:`spatialtis.plotting.base.dot_matrix`,
            :class:`spatialtis.plotting.base.heatmap`,
            :class:`spatialtis.plotting.base.graph_layout_static` or
            :class:`spatialtis.plotting.base.graph_layout_interactive`

    {pyecharts_tips}

    """
    if key is None:
        key = ANALYSIS["spatial_enrichment_analysis"].last_used_key

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]
    order = params["order"]

    if groupby is None:
        groupby = [exp_obs[0]]

    if selected_markers is not None:
        df = df[
            df["marker1"].isin(selected_markers) & df["marker2"].isin(selected_markers)
        ]
    unique_types = np.unique(df["marker1"])
    if use == "heatmap":
        df = df.pivot(index=exp_obs, columns=["marker1", "marker2"], values="value")
        old_cols = df.columns.tolist()
        new_cols = []
        # resort the data
        if selected_markers is None:
            selected_markers = unique_types
        if order:
            for i in product(selected_markers, repeat=2):
                if i in old_cols:
                    new_cols.append(i)
        else:
            for i in combinations_with_replacement(selected_markers, 2):
                if i in old_cols:
                    new_cols.append(i)
        sort_index = pd.MultiIndex.from_tuples(new_cols, names=["marker1", "marker2"])
        df = df.loc[:, sort_index]
        plot_kwargs = dict(
            categorical_colorbar=["Avoidance", "Association"],
            palette=["#2f71ab", "#f7f7f7", "#ba262b"],
        )
        unique_values = np.unique(df.to_numpy())
        if 1 not in unique_values:
            plot_kwargs["palette"] = ["#2f71ab", "#f7f7f7"]
            plot_kwargs["categorical_colorbar"] = ["Avoidance", "Non-sign"]
        if -1 not in unique_values:
            plot_kwargs["palette"] = ["#f7f7f7", "#ba262b"]
            plot_kwargs["categorical_colorbar"] = ["Non-sign", "Association"]
        for k, v in kwargs.items():
            plot_kwargs[k] = v
        return heatmap(
            df,
            row_colors=groupby,
            col_colors=["marker1", "marker2"],
            clustermap_kwargs=dict(row_cluster=True, col_cluster=False),
            **plot_kwargs,
            saved_name="spatial_enrichment",
        )
    else:

        def count(v):
            return {0: 0, 1: 0, -1: 0, **Counter(v)}

        ndf = df.pivot_table(
            index=["marker1", "marker2"], values="value", aggfunc=count
        )
        plot_data = []
        for ix, r in ndf.iterrows():
            r = r["value"]
            no, association, avoidance = r[0], r[1], r[-1]
            sum_all = no + association + avoidance
            sign_size = (association + avoidance) / sum_all
            sign_dir = (
                0
                if sign_size == 0
                else (association if association > avoidance else -avoidance)
                / (association + avoidance)
            )
            plot_data.append([sign_size, sign_dir])
        plot_df = pd.DataFrame(
            data=plot_data, columns=["size", "color"], index=ndf.index
        ).reset_index()

        if use == "dot":
            plot_size = plot_df.pivot(index="marker1", columns="marker2", values="size")
            plot_color = plot_df.pivot(
                index="marker1", columns="marker2", values="color"
            )

            order = plot_size.count().sort_values().index.tolist()
            plot_size = plot_size.loc[order, order[::-1]]
            plot_color = plot_color.loc[order, order[::-1]]

            plot_kwargs = dict(
                show_ticks=False,
                size_legend_title="Sign' ROI",
                matrix_cbar_text=["-1", "1"],
                dot_cbar_text=["Avoidance", "Association"],
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v

            return dot_matrix(
                dot_size=plot_size,
                dot_color=plot_color,
                xlabels=plot_size.columns,
                ylabels=plot_color.index,
                **plot_kwargs,
                saved_name="spatial_enrichment",
            )
        elif use == "graph_interactive":
            plot_kwargs = dict(
                max_width=10,
                width_vmin=0,
                width_vmax=1,
                color_vmin=-1,
                color_vmax=1,
                directed=order,
                cbar_text=["Association", "Avoidance"],
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            return graph_layout_interactive(
                plot_df.to_numpy().tolist(),
                **plot_kwargs,
                saved_name="spatial_enrichment",
            )
        else:
            plot_kwargs = dict(
                max_width=10,
                width_vmin=0,
                width_vmax=1,
                color_vmin=-1,
                color_vmax=1,
                directed=order,
                cbar_text=["Avoidance", "Association"],
                layout="circular_layout",
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            return graph_layout_static(
                plot_df.to_numpy().tolist(),
                **plot_kwargs,
                saved_name="spatial_enrichment",
            )
