from typing import Optional

import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG

from ..utils import timer
from ._neighbors import Neighbors
from ._util import check_neighbors


@timer(prefix="Running community detection")
def communities(
    n: Neighbors,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """leidenalg algorithm for communities detection

    Args:
        n: plotting.Neighbors instance
        export: whether to export the result to anndata.obs
        export_key: the key used to export
        return_df: whether to return the result

    """
    try:
        import leidenalg
    except ImportError:
        raise ImportError("Required leidenalg, try pip install leidenalg.")

    if export_key is None:
        export_key = CONFIG.community_key
    else:
        CONFIG.community_key = export_key

    check_neighbors(n)

    sub_comm = []
    graphs = n.to_graphs()
    for _, graph in tqdm(graphs.items(), **CONFIG.tqdm(desc="find communities"),):
        part = leidenalg.find_partition(
            graph, leidenalg.CPMVertexPartition, resolution_parameter=0.05
        )
        sub_comm += part.membership

    sub_comm = pd.Series(sub_comm, index=n.data.index)
    if export:
        n.adata.obs[export_key] = sub_comm

    if return_df:
        return sub_comm
