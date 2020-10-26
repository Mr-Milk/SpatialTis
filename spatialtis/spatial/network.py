from typing import Any, Optional

import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import check_neighbors
from spatialtis.utils import get_default_params, reuse_docstring, timer


@timer(prefix="Running community detection")
@get_default_params
@reuse_docstring()
def communities(
    n: Neighbors,
    partition_type: Optional[Any] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    **kwargs,
):
    """Leidenalg algorithm for communities detection

    Args:
        n: {n}
        partition_type: The leidenalg partition type
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        **kwargs: Pass to leidenalg.find_partition

    """
    try:
        import leidenalg
    except ImportError:
        raise ImportError("Required leidenalg, try pip install leidenalg.")

    if partition_type is None:
        partition_type = leidenalg.CPMVertexPartition

    if export_key is None:
        export_key = CONFIG.community_key
    else:
        CONFIG.community_key = export_key

    check_neighbors(n)

    partition_kwargs = dict(resolution_parameter=0.05)
    for k, v in kwargs.items():
        partition_kwargs[k] = v

    sub_comm = []
    graphs = n.to_graphs()
    for _, graph in tqdm(graphs.items(), **CONFIG.tqdm(desc="Communities detection"),):
        part = leidenalg.find_partition(graph, partition_type, **partition_kwargs)
        sub_comm += part.membership

    sub_comm = pd.Series(sub_comm, index=n.data.index)
    if export:
        n.adata.obs[export_key] = sub_comm

    if return_df:
        return sub_comm
