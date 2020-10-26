from typing import Optional, Sequence

import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import get_colors, sankey
from spatialtis.utils import adata_uns2df, reuse_docstring


@reuse_docstring()
def exp_neighcells(
    adata: AnnData,
    key: Optional[str] = None,
    score: float = 0.5,
    palette: Optional[Sequence] = None,
    **kwargs,
):
    """(pyecharts) plotting function for expression influenced by neighbor cells

    Args:
        adata: {adata_plotting}
        key: {key}
        score: Threshold for score
        palette: Control the color
        **kwargs: Pass to `Sankey plot <plotting.html#spatialtis.plotting.sankey>`_

    """
    if key is None:
        key = CONFIG.exp_neighcell_key

    df = adata_uns2df(adata, key)
    df = df[df["Score"] >= score]

    cell_types = pd.unique(df.iloc[:, [0, 2]].values.flatten())
    gene_types = pd.unique(df.iloc[:, [1]].values.flatten())

    if palette is None:
        palette = ["Set3", "Spectral"]

    colors = get_colors(len(cell_types) + len(gene_types), palette)
    cell_colormap = dict(zip(cell_types, colors[0 : len(cell_types)]))
    gene_colormap = dict(
        zip(gene_types, colors[len(cell_types) : (len(cell_types) + len(gene_types))])
    )

    nodes = []
    nodes_colors = []
    links = []
    for c in pd.unique(df.iloc[:, [0]].values.flatten()):
        nodes.append(c)
        nodes_colors.append(cell_colormap[c])
    for c in pd.unique(df.iloc[:, [2]].values.flatten()):
        nodes.append(c + " ")
        nodes_colors.append(cell_colormap[c])
    for g in gene_types:
        nodes.append(g)
        nodes_colors.append(gene_colormap[g])
    for i, row in df.iterrows():
        links.append((row[0], row[1], 1))
        links.append((row[1], row[2] + " ", 1))
    p = sankey(nodes, nodes_colors, links, **kwargs)

    # return nodes, nodes_colors, links
    return p
