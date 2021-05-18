from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import pearsonr

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import (
    dot_matrix,
    graph_layout_interactive,
    graph_layout_static,
    heatmap,
)
from spatialtis.utils import doc, get_result


@doc
def neighborhood_analysis(
        data: AnnData,
        groupby: Optional[Sequence[str]] = None,
        selected_types: Optional[Sequence] = None,
        key: Optional[str] = None,
        use: str = "dot_matrix",  # heatmap, dot_matrix
        **kwargs,
):
    """Visualization for neighborhood analysis

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        use: "dot_matrix", "heatmap", "graph_static" or "graph_interactive" (Default: "dot_matrix")
        **kwargs: Pass to :class:`spatialtis.plotting.base.dot_matrix`,
            :class:`spatialtis.plotting.base.heatmap`,
            :class:`spatialtis.plotting.base.graph_layout_static` or
            :class:`spatialtis.plotting.base.graph_layout_interactive`

    {pyecharts_tips}

    """
    if key is None:
        key = ANALYSIS["neighborhood_analysis"].last_used_key

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]
    order = params["order"]

    if groupby is None:
        groupby = [exp_obs[0]]

    if selected_types is not None:
        df = df[df["type1"].isin(selected_types) & df["type2"].isin(selected_types)]
    unique_types = np.unique(df["type1"])
    if use == "heatmap":
        df = df.pivot(index=exp_obs, columns=["type1", "type2"], values="value")
        # resort the data
        if selected_types is None:
            selected_types = unique_types
        if order:
            cols = [i for i in product(selected_types, repeat=2)]
        else:
            cols = [i for i in combinations_with_replacement(selected_types, 2)]

        old_cols = df.columns.tolist()
        # filp the order in old index
        new_cols = []
        for i in old_cols:
            if i in cols:
                new_cols.append(i)
            else:
                new_cols.append(i[::-1])
        # replace old with new
        df.columns = pd.MultiIndex.from_tuples(new_cols, names=["type1", "type2"])
        sortIndex = dict(zip(cols, range(len(cols))))
        new_cols = sorted(new_cols, key=lambda x: sortIndex[x])
        # sort it
        sort_index = pd.MultiIndex.from_tuples(new_cols, names=["type1", "type2"])
        df = df.loc[:, sort_index]
        plot_kwargs = dict(
            categorical_colorbar=["Avoidance", "Association"],
            palette=["#2f71ab", "#f7f7f7", "#ba262b"],
            clustermap_kwargs=dict(row_cluster=True, col_cluster=False, method="ward"),
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
            col_colors=["type1", "type2"],
            **plot_kwargs,
            saved_name="neighborhood_analysis",
        )
    else:

        def count(v):
            return {0: 0, 1: 0, -1: 0, **Counter(v)}

        ndf = df.pivot_table(index=["type1", "type2"], values="value", aggfunc=count)
        plot_data = []
        for ix, r in ndf.iterrows():
            r = r["value"]
            no, association, avoidance = r[0], r[1], r[-1]
            sum_all = no + association + avoidance
            real_size = association + avoidance
            sign_size = real_size / sum_all
            sign_dir = (
                0
                if sign_size == 0
                else (association if association > avoidance else -avoidance)
                     / (association + avoidance)
            )
            plot_data.append([real_size, sign_size, sign_dir])
        plot_df = pd.DataFrame(
            data=plot_data, columns=["real_size", "size", "color"], index=ndf.index
        ).reset_index()
        if use == "dot_matrix":
            try:
                df = get_result(data, ANALYSIS["cell_components"].last_used_key)
                cc = {}
                for n, g in df.groupby("type"):
                    cc[n] = g["value"].tolist()
            except KeyError:
                raise Exception(
                    "Please run cell_components before plotting the dot matrix plot"
                )
            corr = {}
            for i in product(unique_types, repeat=2):
                corr[i] = [pearsonr(cc[i[0]], cc[i[1]])[0]]

            ccdf = pd.DataFrame(corr).T.reset_index()
            ccdf.columns = ["type1", "type2", "corr"]
            plot_df = (
                plot_df.set_index(["type1", "type2"])
                    .join(
                    ccdf.set_index(["type1", "type2"]),
                    on=["type1", "type2"],
                    how="left",
                )
                    .reset_index()
            )

            plot_corr = plot_df.pivot(index="type1", columns="type2", values="corr")
            plot_size = plot_df.pivot(
                index="type1", columns="type2", values="real_size"
            )
            plot_color = plot_df.pivot(index="type1", columns="type2", values="color")

            order = plot_corr.count().sort_values().index.tolist()
            plot_corr = plot_corr.loc[order, order[::-1]]
            plot_size = plot_size.loc[order, order[::-1]]
            plot_color = plot_color.loc[order, order[::-1]]

            plot_kwargs = dict(
                show_ticks=False,
                size_legend_title="Sign' ROI",
                matrix_cbar_title="■\nPearson\nCorrelation",
                matrix_cbar_text=["-1", "1"],
                dot_cbar_text=["Avoidance", "Association"],
                dot_cbar_title="●\n% of\ninteraction",
                xtickslabel_rotation=90,
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v

            return dot_matrix(
                matrix=plot_corr,
                dot_size=plot_size,
                dot_color=plot_color,
                xlabels=plot_corr.columns,
                ylabels=plot_corr.index,
                **plot_kwargs,
                saved_name="neighborhood_analysis",
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
                saved_name="neighborhood_analysis",
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
                saved_name="neighborhood_analysis",
            )
