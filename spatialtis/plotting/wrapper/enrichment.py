from collections import Counter
from itertools import product
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import norm

from spatialtis.config import CONFIG
from spatialtis.plotting.base import DotMatrix, TriDotMatrix, heatmap
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def spatial_enrichment_analysis(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    use: str = "dot_matrix",
    key: Optional[str] = None,
    **kwargs,
):
    """(matplotlib) plotting function for plotting enrichment analysis

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        use: Options are "dot_matrix" and "heatmap"
        key: {key}
        **kwargs: Pass to `Dot-matrix plot <plotting.html#spatialtis.plotting.dot_matrix>`_ or
                  `Heatmap <plotting.html#spatialtis.plotting.heatmap>`_

    """
    if key is None:
        key = CONFIG.spatial_enrichment_analysis_key

    df, params = adata_uns2df(adata, key, params=True)
    order = True  # params["order"]
    pval = params["pval"]

    data = df.to_numpy()
    sign_data = []
    for arr in data:
        sign_arr = []
        for i in arr:
            p_value = norm.sf(abs(i)) * 2
            if p_value >= pval:
                sign_arr.append(0)
            else:
                if i > 0:
                    sign_arr.append(1)
                else:
                    sign_arr.append(-1)
        sign_data.append(sign_arr)
    df = pd.DataFrame(sign_data, columns=df.columns, index=df.index)

    if selected_types is None:
        selected_types = pd.unique(df.columns.to_frame(index=False).iloc[:, 0])
    combs = [i for i in product(selected_types, repeat=2)]
    cols = df.columns.tolist()
    colnames = df.columns.names
    new_cols = []
    for comb in combs:
        if comb not in cols:
            if comb[::-1] in cols:
                new_cols.append(comb[::-1])
        else:
            new_cols.append(comb)
    df = df[new_cols]
    df.columns = pd.MultiIndex.from_tuples(combs)
    df.columns.names = colnames

    if use == "dot_matrix":

        # [cell 1, cell 2, interaction type, No. of sig (both type), % of sign type, p)
        counts = list()
        for label, data in df.items():
            # count the No. of sign for both type
            c_data = sorted(
                [(k, v) for k, v in {**{0: 0, 1: 0, -1: 0}, **Counter(data)}.items()],
                key=lambda x: x[1],
                reverse=True,
            )

            all_sign = 0
            for (t, n) in c_data:
                if t != 0:
                    all_sign += n

            if all_sign == 0:
                counts.append([label[0], label[1], 0, 0])
            else:
                # if most ROI are no relationship
                if c_data[0][0] == 0:
                    # check the next
                    if c_data[1][0] == 1:
                        counts.append(
                            [label[0], label[1], all_sign, c_data[1][1] / all_sign]
                        )
                    else:
                        counts.append(
                            [label[0], label[1], all_sign, 1 - c_data[1][1] / all_sign]
                        )
                # if most ROI are interaction
                elif c_data[0][0] == 1:
                    counts.append(
                        [label[0], label[1], all_sign, c_data[0][1] / all_sign]
                    )
                else:
                    counts.append(
                        [label[0], label[1], all_sign, 1 - c_data[0][1] / all_sign]
                    )
        counts = pd.DataFrame(counts, columns=["type1", "type2", "all", "%"])
        dot_color = counts.pivot(index="type1", columns="type2", values="%")
        dot_color = dot_color.loc[selected_types, selected_types]
        dot_size = counts.pivot(index="type1", columns="type2", values="all")
        dot_size = dot_size.loc[selected_types, selected_types]
        xlabels = dot_size.columns
        ylabels = dot_size.index
        dot_color = dot_color.to_numpy()
        dot_size = dot_size.to_numpy()

        def filter_nan(matrix):
            diag_matrix = []
            for i in matrix:
                arr = []
                for t in i:
                    if not np.isnan(t):
                        arr.append(t)
                diag_matrix.append(arr)
            return diag_matrix

        if order:
            plot_kwargs = dict(
                show_ticks=False,
                size_legend_title="Sign' ROI",
                dot_cbar_text=["Avoidance", "Association"],
                dot_cbar_title="●\n% of\ninteraction",
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            p = DotMatrix(
                dot_color=dot_color,
                dot_size=dot_size,
                xlabels=xlabels,
                ylabels=ylabels,
                **plot_kwargs,
            )
        else:
            diag_color = filter_nan(dot_color)
            diag_size = filter_nan(dot_size)

            plot_kwargs = dict(
                show_ticks=False,
                size_legend_title="Sign' ROI",
                dot_cbar_text=["Avoidance", "Association"],
                dot_cbar_title="●\n% of\ninteraction",
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            p = TriDotMatrix(
                diagonal_dot_color=diag_color,
                diagonal_dot_size=diag_size,
                labels=xlabels,
                **plot_kwargs,
            )

    else:
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors=["Marker1", "Marker2"],
            palette=["#2f71ab", "#f7f7f7", "#ba262b"],
            colorbar_type="categorical",
            categorical_colorbar_text=["No co-expression", "Co-expression"],
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
