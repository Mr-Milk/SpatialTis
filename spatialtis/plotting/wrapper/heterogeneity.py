from typing import Optional, Sequence

from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import stacked_bar, violin_plot
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def spatial_heterogeneity(
    adata: AnnData,
    groupby: Optional[Sequence[str]] = None,
    key: Optional[str] = None,
    metric: str = "heterogeneity",
    **kwargs,
):
    """(bokeh) plotting function for cell morphology

    Args:
        adata: {adata_plotting}
        groupby: {groupby}
        key: {key}
        metric: "heterogeneity" or "KL", "KL" only available if you use shannon entropy
        **kwargs: Pass to `Violin plot <plotting.html#spatialtis.plotting.violin_plot>`_ or
                  `Stacked-bar plot <plotting.html#spatialtis.plotting.stacked_bar>`_

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
