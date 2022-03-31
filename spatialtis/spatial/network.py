from typing import Any, Dict, Optional

import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import euclidean

from spatialtis.abc import AnalysisBase
from spatialtis.utils import col2adata_obs, doc
from spatialtis.utils import try_import, read_neighbors


@doc
def cell_community(data: AnnData,
                   resolution: float = 0.05,
                   partition_type: Optional[Any] = None,
                   partition_kwargs: Optional[Dict] = None,
                    export_key: str = "community_id",
                   **kwargs, ):
    """Spatial communities detection

    Here we use Leiden graph cluster algorithm

    Args:
        data: {adata}
        resolution:
        partition_type: The leidenalg partition type
        partition_kwargs: Pass to leidenalg.find_partition
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    """

    ab = AnalysisBase(data,
                      display_name="Cell community",
                      export_key=export_key,
                      **kwargs)

    # import leidenalg
    # import igraph as ig
    leidenalg = try_import("leidenalg")
    ig = try_import("igraph", install_name="python-igraph")

    ab.check_neighbors()

    if partition_type is None:
        partition_type = leidenalg.CPMVertexPartition
    if partition_kwargs is None:
        partition_kwargs = {"resolution_parameter": resolution}
    else:
        partition_kwargs = {"resolution_parameter": 0.05, **partition_kwargs}

    graphs = []
    track_ix = []
    sub_comm = []
    for roi_name, roi_data, points in ab.roi_iter_with_points():
        labels = roi_data[ab.cell_id_key]
        neighbors = read_neighbors(roi_data, ab.neighbors_key)
        vertices = []
        edge_mapper = {}
        for i, (x, y) in zip(labels, points):
            vertices.append({"name": i, "x": x, "y": y})
            edge_mapper[i] = (x, y)

        graph_edges = []
        for k, vs in zip(labels, neighbors):
            if len(vs) > 0:
                for v in vs:
                    if k < v:
                        distance = euclidean(edge_mapper[k], edge_mapper[v])
                        graph_edges.append(
                            {"source": k, "target": v, "weight": distance}
                        )
        graph = ig.Graph.DictList(vertices, graph_edges)
        part = leidenalg.find_partition(graph, partition_type, **partition_kwargs)
        sub_comm += part.membership
        graphs.append(graph)
        track_ix.append(roi_name)

    sub_comm = pd.Series(sub_comm, index=data.obs.index)
    col2adata_obs(sub_comm, data, ab.export_key)
    ab.stop_timer()
