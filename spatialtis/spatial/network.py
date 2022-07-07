import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import euclidean
from typing import Any, Dict

from spatialtis.abc import AnalysisBase
from spatialtis.utils import col2adata, doc
from spatialtis.utils import try_import


@doc
def cell_community(data: AnnData,
                   resolution: float = 0.05,
                   partition_type: Any = None,
                   partition_kwargs: Dict = None,
                   export_key: str = "community_id",
                   **kwargs, ):
    """Spatial communities detection

    Here we use Leiden graph cluster algorithm

    Parameters
    ----------
    data : {adata}
    resolution : float, default: 0.05
        Control the process of partition.
    partition_type :
        The leidenalg partition type.
    partition_kwargs :
        Pass to leidenalg.find_partition.
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

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
    for roi_name, labels, neighbors, points in ab.iter_roi(fields=['neighbors', 'centroid']):
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
    col2adata(sub_comm, data, ab.export_key)
    ab.stop_timer()
