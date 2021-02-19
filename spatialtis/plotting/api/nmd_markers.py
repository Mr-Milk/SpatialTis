from typing import Optional

from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import graph_layout_static
from spatialtis.utils import doc, get_result


@doc
def NMDMarkers(
    adata: AnnData, key: Optional[str] = None, **kwargs,
):
    """(pyecharts) plotting function for expression influenced by neighbor cells

    Args:
        adata: {adata_plotting}
        key: {key}
        **kwargs: Pass to :class:`spatialtis.plotting.base.graph_layout_static`

    """
    if key is None:
        key = ANALYSIS["NMDMarkers"].last_used_key

    df = get_result(adata, key)
    df = df[["neighbor_marker", "marker", "dependency", "corr"]]
    plot_kwargs = dict(
        width_vmin=0,
        width_vmax=1,
        color_vmin=-1,
        color_vmax=1,
        layout="circular_layout",
        directed=True,
    )
    for k, v in kwargs.items():
        plot_kwargs[k] = v
    return graph_layout_static(df.to_numpy(), **plot_kwargs, saved_name="NMDMarkers",)
