from collections import Counter
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import dotplot, heatmap
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def spatial_distribution(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    selected_types: Optional[Sequence] = None,
    key: Optional[str] = None,
    use: str = "dot",  # heatmap
    **kwargs,
):
    """(matplotlib) plotting function for plotting distribution

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        selected_types: {selected_types}
        key: {key}
        use: Options are "dot" and "heatmap"
        **kwargs: Pass to `Dotplot <plotting.html#spatialtis.plotting.dotplot>`_ or
                  `Heatmap <plotting.html#spatialtis.plotting.heatmap>`_

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
        tb = tb.drop(tb.tail(1).index).fillna(0)
        tb.index = names
        tb = tb[[0, 1, 2, 3]]
        colors = np.array(["#FFC408", "#c54a52", "#4a89b9", "#5a539d"] * len(tb))

        plot_kwargs = dict(legend_title="ROI", annotated=False, alpha=1,)

        for k, v in kwargs.items():
            plot_kwargs[k] = v

        p = dotplot(
            tb.to_numpy(),
            colors,
            xlabels=["No Cell", "Random", "Regular", "Cluster"],
            ylabels=names,
            **plot_kwargs,
        )
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
