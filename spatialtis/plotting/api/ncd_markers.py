from typing import Optional

from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import graph_layout_interactive, graph_layout_static
from spatialtis.utils import doc, get_result


@doc
def NCDMarkers(
        adata: AnnData, key: Optional[str] = None, use: str = "static", **kwargs,
):
    """Plotting function for expression influenced by neighbor cells

    Args:
        adata: {adata_plotting}
        key: {key}
        use: "static" or "interactive" (Default: "static")
        **kwargs: Pass to :class:`spatialtis.plotting.base.graph_layout_static`

    """
    if key is None:
        key = ANALYSIS["NCDMarkers"].last_used_key

    df = get_result(adata, key)
    if "cell_type" in df.columns:
        df["marker"] = df["cell_type"] + df["marker"]
    df["color"] = 1
    df = df[["neighbor_type", "marker", "dependency", "color"]]
    plot_kwargs = dict(
        width_vmin=0,
        width_vmax=1,
        color_vmin=-1,
        color_vmax=1,
        layout="bipartite_layout",
        curve=0,
        directed=True,
    )
    for k, v in kwargs.items():
        plot_kwargs[k] = v
    if use == "interactive":
        return graph_layout_interactive(
            df.to_numpy(), **plot_kwargs, saved_name="NCDMarkers",
        )
    else:
        return graph_layout_static(
            df.to_numpy(), **plot_kwargs, saved_name="NCDMarkers",
        )
