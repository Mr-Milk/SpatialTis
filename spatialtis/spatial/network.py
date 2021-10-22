from ast import literal_eval
from typing import Any, Dict, Optional

import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import euclidean

from spatialtis.abc import AnalysisBase
from spatialtis.config import Config
from spatialtis.utils import NeighborsNotFoundError
from spatialtis.utils import col2adata_obs, doc, pbar_iter


@doc
class cell_community(AnalysisBase):
    """Spatial communities detection

    Here we use Leiden graph cluster algorithm

    Args:
        data: {adata}
        partition_type: The leidenalg partition type
        partition_kwargs: Pass to leidenalg.find_partition
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        partition_type: Optional[Any] = None,
        partition_kwargs: Optional[Dict] = None,
        **kwargs,
    ):
        super().__init__(data, **kwargs)

        try:
            import leidenalg
        except ImportError:
            raise ImportError("Required leidenalg, try pip install leidenalg.")

        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "Required python-igraph, try `pip install python-igraph`."
            )

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if partition_type is None:
            partition_type = leidenalg.CPMVertexPartition
        if partition_kwargs is None:
            partition_kwargs = {"resolution_parameter": 0.05}
        else:
            partition_kwargs = {"resolution_parameter": 0.05, **partition_kwargs}

        need_eval_cent = self.is_col_str(self.centroid_key)
        need_eval_neigh = self.is_col_str(self.neighbors_key)

        graphs = []
        names = []
        for n, g in data.obs.groupby(self.exp_obs):
            if need_eval_cent:
                centroids = [literal_eval(c) for c in g[self.centroid_key]]
            else:
                centroids = [c for c in g[self.centroid_key]]
            if need_eval_neigh:
                neighbors = [literal_eval(n) for n in g[self.neighbors_key]]
            else:
                neighbors = [n for n in g[self.neighbors_key]]

            vertices = []
            edge_mapper = {}
            for i, (x, y) in zip(g[Config.neighbors_ix_key], centroids):
                vertices.append({"name": i, "x": x, "y": y})
                edge_mapper[i] = (x, y)

            graph_edges = []
            for k, vs in zip(g[Config.neighbors_ix_key], neighbors):
                if len(vs) > 0:
                    for v in vs:
                        if k != v:
                            distance = euclidean(edge_mapper[k], edge_mapper[v])
                            graph_edges.append(
                                {"source": k, "target": v, "weight": distance}
                            )
            graphs.append(ig.Graph.DictList(vertices, graph_edges))
            names.append(n)

        neighbors_graphs = dict(zip(names, graphs))
        sub_comm = []
        for _, graph in pbar_iter(
            neighbors_graphs.items(),
            desc="Communities detection",
        ):
            part = leidenalg.find_partition(graph, partition_type, **partition_kwargs)
            sub_comm += part.membership

        sub_comm = pd.Series(sub_comm, index=data.obs.index)
        col2adata_obs(sub_comm, self.data, self.export_key)
        self.stop_timer()
