from collections import Counter
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import dot_plot, heatmap
from spatialtis.utils import doc, get_result


@doc
def spatial_distribution(
    data: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "dot",  # heatmap
    **kwargs,
):
    """Visualization for spatial distribution

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        use: "dot" or "heatmap" (Default: "dot")
        **kwargs: Pass to :class:`spatialtis.plotting.base.dot_matrix` or
            :class:`spatialtis.plotting.base.heatmap`

    """
    if key is None:
        key = ANALYSIS["spatial_distribution"].last_used_key

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]

    if groupby is None:
        groupby = exp_obs

    if selected_types is not None:
        df = df[df["type"].isin(selected_types)]

    if use == "dot":

        plot_df = pd.DataFrame(
            [
                {0: 0, 1: 0, 2: 0, 3: 0, **Counter(g["pattern"])}
                for _, g in df.groupby("type")
            ]
        ).rename(
            columns=dict(
                zip([0, 1, 2, 3], ["No cells", "Random", "Regular", "Cluster"])
            )
        )

        colors = np.array(["#FFC408", "#c54a52", "#4a89b9", "#5a539d"] * len(plot_df))

        plot_kwargs = dict(legend_title="ROI", alpha=1, xtickslabel_rotation=90)
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        return dot_plot(
            plot_df, colors, **plot_kwargs, saved_name="spatial_distribution",
        )

    else:
        plot_df = df.pivot(index=exp_obs, columns="type", values="pattern").fillna(0)
        plot_kwargs = dict(
            row_colors=groupby,
            col_colors="type",
            palette=["#fffec6", "#c54a52", "#4a89b9", "#5a539d"],
            categorical_colorbar=["No Cell", "Random", "Regular", "Cluster"],
            clustermap_kwargs=dict(row_cluster=None, col_cluster=True,),
        )
        # allow user to overwrite the default plot config
        for k, v in kwargs.items():
            plot_kwargs[k] = v

        return heatmap(plot_df, **plot_kwargs, saved_name="spatial_distribution",)
