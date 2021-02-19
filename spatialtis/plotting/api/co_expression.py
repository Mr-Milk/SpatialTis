from typing import Optional, Sequence

import pandas as pd
from anndata import AnnData

from spatialtis import get_result
from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import (
    graph_layout_interactive,
    graph_layout_static,
    heatmap,
)
from spatialtis.utils import doc


@doc
def spatial_co_expression(
    data: AnnData,
    selected_markers: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "graph_static",  # heatmap
    **kwargs,
):
    """Visualization for spatial enrichment analysis

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_markers: {selected_markers}
        key: {key}
        use: "heatmap", "graph_static" or "graph_interactive" (Default: "graph_static")
        **kwargs: Pass to :class:`spatialtis.plotting.base.heatmap`,
            :class:`spatialtis.plotting.base.graph_layout_static` or
            :class:`spatialtis.plotting.base.graph_layout_interactive`

    {pyecharts_tips}

    """
    if key is None:
        key = ANALYSIS["spatial_co_expression"].last_used_key

    df = get_result(data, key)
    if selected_markers is not None:
        df = df[
            (df["marker1"].isin(selected_markers))
            & (df["marker2"].isin(selected_markers))
        ]

    if use == "heatmap":
        ndf = df.copy()
        ndf.columns = ["marker2", "marker1", "corr", "pvalue"]
        plot_df = (
            pd.concat([df, ndf])
            .drop_duplicates()
            .pivot(columns="marker1", index="marker2", values="corr")
            .fillna(0)
        )
        return heatmap(plot_df, row_label="marker2", col_label="marker1", **kwargs)
    else:
        df["width"] = 1
        df = df[["marker1", "marker2", "width", "corr"]]
        plot_kwargs = dict(
            width_vmin=0, color_vmin=-1, color_vmax=1, saved_name="co-expression"
        )
        for k, v in kwargs.items():
            plot_kwargs[k] = v
        if use == "graph_interactive":
            return graph_layout_interactive(
                df.to_numpy(), cbar_text=["pos", "neg"], **plot_kwargs
            )
        else:
            return graph_layout_static(
                df.to_numpy(), cbar_text=["neg", "pos"], **plot_kwargs
            )
