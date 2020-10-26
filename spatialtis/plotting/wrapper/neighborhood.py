from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import norm, pearsonr

from spatialtis.config import CONFIG
from spatialtis.plotting.base import DotMatrix, TriDotMatrix, heatmap
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def neighborhood_analysis(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "dot_matrix",  # heatmap, dot_matrix
    **kwargs,
):
    """("dot_matrix", "heatmap": matplotlib) plotting function for neighborhood analysis

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        use: Options are "dot_matrix" and "heatmap"
        **kwargs: Pass to `Dot-matrix plot <plotting.html#spatialtis.plotting.dot_matrix>`_ or
                  `Heatmap <plotting.html#spatialtis.plotting.heatmap>`_

    """
    if key is None:
        key = CONFIG.neighborhood_analysis_key

    df, params = adata_uns2df(adata, key, params=True)
    order = params["order"]
    method = params["method"]
    pval = params["pval"]

    if method == "zscore":
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
        selected_types = np.unique(df.columns.to_frame(index=False).to_numpy())
    if order:
        combs = [i for i in product(selected_types, repeat=2)]
        df = df[combs]
    else:
        combs = [i for i in combinations_with_replacement(selected_types, 2)]
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

    """
    if use == "graph":
        p = cc_interactions(df, {-1: "Avoidance", 1: "Association"}, **kwargs)
    """
    if use == "dot_matrix":
        try:
            cc = adata_uns2df(adata, CONFIG.cell_components_key)
        except KeyError:
            raise Exception(
                "Please run cell_components before plotting the dot matrix plot"
            )

        # [cell 1, cell 2, interaction type, No. of sig (both type), % of sign type, p)
        counts = list()
        for label, data in df.items():
            # count pearson correlation between cell components
            p = pearsonr(cc[label[0]], cc[label[1]])[0]

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
                counts.append([label[0], label[1], 0, 0, p])
            else:
                # if most ROI are no relationship
                if c_data[0][0] == 0:
                    # check the next
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
                # if most ROI are interaction
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
        # the pivot operation will sort alphabetically, force to sort again
        matrix = matrix.loc[selected_types, selected_types]
        xlabels = matrix.columns
        ylabels = matrix.index
        matrix = matrix.to_numpy()
        dot_color = counts.pivot(index="type1", columns="type2", values="%")
        dot_color = dot_color.loc[selected_types, selected_types].to_numpy()
        dot_size = counts.pivot(index="type1", columns="type2", values="all")
        dot_size = dot_size.loc[selected_types, selected_types].to_numpy()

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
                matrix_cbar_title="■\nPearson\nCorrelation",
                matrix_cbar_text=["-1", "1"],
                dot_cbar_text=["Avoidance", "Association"],
                dot_cbar_title="●\n% of\ninteraction",
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            p = DotMatrix(
                matrix,
                dot_color,
                dot_size,
                xlabels=xlabels,
                ylabels=ylabels,
                **plot_kwargs,
            )
        else:
            diag_matrix = filter_nan(matrix)
            diag_color = filter_nan(dot_color)
            diag_size = filter_nan(dot_size)

            plot_kwargs = dict(
                show_ticks=False,
                size_legend_title="Sign' ROI",
                block_cbar_title="■\nPearson\nCorrelation",
                block_cbar_text=["-1", "1"],
                dot_cbar_text=["Avoidance", "Association"],
                dot_cbar_title="●\n% of\ninteraction",
            )
            for k, v in kwargs.items():
                plot_kwargs[k] = v
            p = TriDotMatrix(
                diag_matrix, diag_color, diag_size, labels=xlabels, **plot_kwargs
            )

    else:
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
