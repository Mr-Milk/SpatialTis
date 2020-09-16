import warnings
from typing import Optional, Sequence, Union

import numpy as np
from anndata import AnnData
from neighborhood_analysis import get_bbox, get_bbox_neighbors, get_point_neighbors
from scipy.spatial.distance import euclidean
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import col2adata_obs, lprint, timer


class Neighbors(object):
    """Storage object for cell relationship

    Args:
        adata: anndata object to perform analysis
        geom: how to resolve cell data, as "point" or "shape"
        groupby: how your experiments grouped, (Default: spatialtis.CONFIG.EXP_OBS)
        type_key: the key name of cell type in anndata.obs (Default: spatialtis.CONFIG.CELL_TYPE_KEY)
        shape_key: the key name of cell shape in anndata.obs (Default: spatialtis.CONFIG.SHAPE_KEY)
        centroid_key: anndata.obs key that store cell centroid info (Default: spatialtis.CONFIG.CENTROID_KEY)

    """

    @timer(verbose=False)
    def __init__(
        self,
        adata: AnnData,
        geom: str = "shape",
        *,
        groupby: Union[Sequence, str, None] = None,
        type_key: Optional[str] = None,
        shape_key: Optional[str] = None,
        centroid_key: Optional[str] = None,
    ):

        # keys for query info from anndata
        if geom not in ["shape", "point"]:
            raise ValueError("Available options for 'geom' are 'shape', 'point'")

        # verify if neighbors and polycells built
        self.__polycells = False
        self.__neighborsbuilt = False

        self.__adata = adata
        self.__data = adata.obs
        self.__geom = geom
        self.__typekey = type_key
        self.__shapekey = shape_key
        self.__centkey = centroid_key
        self.__groupby = groupby
        self.__groups = self.__data.groupby(groupby, sort=False)

        # define vars for storage
        self.__polycellsdb = None
        self.__neighborsdb = None
        self.__uniquetypes = np.unique(self.__data[type_key])
        self.__types = {n: list(g[self.__typekey]) for n, g in self.__groups}

    def __repr__(self):
        return "A spatialtis.Neighbors instance"

    @timer(prefix="Finding cell neighbors")
    def find_neighbors(
        self, expand: Optional[float] = None, scale: Optional[float] = None,
    ):
        """To find the neighbors of each cell

        Args: expand: If the cell is shape, it means how much units to expand each cell; If the cell is point, it's the search radius
        scale: how much to scale each cell, only if cell is shape

        """
        # handle parameters
        if (expand is None) & (scale is None):
            scale = 1.0
            expand = None
            warnings.warn("Neither 'expand' or 'scale' are specific")
        elif (expand is not None) & (scale is not None):
            warnings.warn(
                f"Conflict parameters, can't set 'expand' and 'scale' in the same time, use expand={expand}"
            )
        elif (expand is None) & (self.__geom == "point"):
            raise ValueError("Parameter 'expand' is not specific")
        elif expand is not None:
            if expand < 0:
                raise ValueError("Can't shrink cell, 'expand' must >= 0")
        elif scale is not None:
            if scale < 1:
                raise ValueError("Can't shrink cell, 'scale' must >= 1")

        if self.__geom == "point":
            lprint("Cell resolved as point data, searching neighbors using KD-tree")
        if self.__geom == "shape":
            lprint("Cell resolved as shape data, searching neighbors using R-tree")

        # prepare data for shape neighbor search
        if (self.__geom == "shape") & (not self.__polycells):
            results = []
            names = []
            for n, g in tqdm(self.__groups, **CONFIG.tqdm(desc="Get cells' bbox",),):
                shapes = g[self.__shapekey]
                polycells = get_bbox([eval(cell) for cell in shapes])
                results.append(polycells)
                names.append(n)

            self.__polycellsdb = dict(zip(names, results))
            self.__polycells = True

        results = []
        names = []
        # shape neighbor search
        if self.__geom == "shape":
            for n, polycells in tqdm(
                self.__polycellsdb.items(), **CONFIG.tqdm(desc="Find neighbors"),
            ):
                if expand is not None:
                    results.append(get_bbox_neighbors(polycells, expand=expand))
                else:
                    results.append(get_bbox_neighbors(polycells, scale=scale))
                names.append(n)
        # point neighbor search
        else:
            for n, g in tqdm(self.__groups, **CONFIG.tqdm(desc="Find neighbors"),):
                cells = [eval(c) for c in g[self.__centkey]]
                results.append(get_point_neighbors(cells, expand))
                names.append(n)

        self.__neighborsdb = dict(zip(names, results))
        self.__neighborsbuilt = True

        return self

    def export_neighbors(self, export_key: Optional[str] = None):
        """Export computed neighbors

        The neighbors relationship are stored in dict, number is index of cell in each ROI
        Correspond to the order in anndata object

        This will export to anndata object's obs field

        Args:
            export_key: the key name to export

        """

        if export_key is None:
            export_key = CONFIG.neighbors_key
        else:
            CONFIG.neighbors_key = export_key

        if not self.__neighborsbuilt:
            return "Please run .find_neighbors() before further analysis."

        neighbors = []
        for n, g in self.__groups:
            neighdict = self.__neighborsdb[n]
            for i in range(0, len(g)):
                try:
                    neighs = neighdict[i]
                except KeyError:
                    neighs = []
                neighbors.append(str(neighs))
        col2adata_obs(neighbors, self.__adata, export_key)

        return self

    def neighbors_count(
        self, export_key: Optional[str] = None,
    ):
        """Get how many neighbors for each cell

        This will write to anndata object's obs field,
        If 'selected_types' are used, this function won't work

        Args:
            export_key: the key name to export

        """

        if export_key is None:
            export_key = CONFIG.neighbors_count_key
        else:
            CONFIG.neighbors_count_key = export_key

        if not self.__neighborsbuilt:
            return "Please run .find_neighbors() before further analysis."

        counts = []
        for n, g in self.__groups:
            neighdict = self.__neighborsdb[n]
            for i in range(0, len(g)):
                try:
                    count = len(neighdict[i])
                except KeyError:
                    count = 0
                counts.append(count)
        col2adata_obs(counts, self.__adata, export_key)

        return self

    def read_neighbors(self, read_key: Optional[str] = None):
        """Read computed neighbors from anndata

        Args:
            read_key: the key name to read

        """

        if read_key is None:
            read_key = CONFIG.neighbors_key

        if read_key not in self.__adata.obs.keys():
            raise KeyError(f"{read_key} not exists.")

        self.__groups = self.__data.groupby(self.__groupby)
        neighborsdb = {}
        for n, g in self.__groups:
            neighborsdb[n] = dict(zip(range(len(g)), [eval(i) for i in g[read_key]]))
        self.__neighborsdb = neighborsdb
        self.__neighborsbuilt = True

        return self

    def to_graphs(self):
        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "Required python-igraph, try `pip install python-igraph`."
            )
        if not self.__neighborsbuilt:
            return None

        graphs = []
        names = []
        for n, g in self.__groups:
            centroids = [eval(c) for c in g[self.__centkey]]
            vertices = [
                {"name": i, "x": x, "y": y} for i, (x, y) in enumerate(centroids)
            ]
            edges = self.__neighborsdb[n]
            graph_edges = []
            for k, vs in edges.items():
                if len(vs) > 0:
                    for v in vs:
                        if k != v:
                            distance = euclidean(centroids[k], centroids[v])
                            graph_edges.append(
                                {"source": k, "target": v, "weight": distance}
                            )
            graphs.append(ig.Graph.DictList(vertices, graph_edges))
            names.append(n)

        new_graphs = dict(zip(names, graphs))

        return new_graphs

    @property
    def neighborsbuilt(self):
        """Check if the neighbors have been built"""
        return self.__neighborsbuilt

    @property
    def unitypes(self):
        """The unique cell types"""
        if self.__typekey is not None:
            return self.__uniquetypes
        else:
            return None

    @property
    def types(self):
        """Cell types order by their index number in anndata"""
        if self.__typekey is not None:
            return self.__types
        else:
            return None

    @property
    def type_key(self):
        """key in anndata.uns use to store cell type info"""
        return self.__typekey

    @property
    def data(self):
        """the info in anndata.obs"""
        return self.__data

    @property
    def adata(self):
        """soft link to anndata object"""
        return self.__adata

    @property
    def neighbors(self):
        """Return the computed neighbors relationship"""
        return self.__neighborsdb

    @property
    def expobs(self):
        """the experiment observations used to process neighbors"""
        return self.__groupby
