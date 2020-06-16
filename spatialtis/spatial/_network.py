from typing import Optional

import leidenalg
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG

from ._neighbors import Neighbors
from ._util import check_neighbors


def communities(
    n: Neighbors, export_key: Optional[str] = None,
):
    """leidenalg algorithm for communities detection

    Args:
        n: spatial.Neighbors instance
        export_key: the key name to store info, exported to anndata.obs field

    """

    if export_key is None:
        export_key = CONFIG.community_key
    else:
        CONFIG.community_key = export_key

    check_neighbors(n)

    sub_comm = []
    graphs = n.to_graphs()
    for _, graph in tqdm(
        graphs.items(),
        desc="find communities",
        bar_format=CONFIG.PBAR_FORMAT,
        disable=(not CONFIG.PROGRESS_BAR),
    ):
        part = leidenalg.find_partition(graph, leidenalg.ModularityVertexPartition)
        sub_comm += part.membership

    sub_comm = pd.Series(sub_comm, index=n.data.index)
    # n.data[export_key] = sub_comm
    n.adata.obs[export_key] = sub_comm
    # n.data[export_key]
