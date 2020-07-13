import warnings
from typing import Optional, Sequence, Union

import igraph as ig
import numpy as np
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.spatial.distance import euclidean
from shapely.affinity import scale as sscale
from shapely.geometry import asMultiPoint, box
from shapely.strtree import STRtree
from shapely.wkt import dumps, loads
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import col2adata_obs, lprint, timer


def _polygonize_cells(shapekey, group):
    shapes = group[shapekey]
    polycells = []
    for cell in shapes:
        c = asMultiPoint(eval(cell))
        polycells.append(c)
    polycells = [dumps(c) for c in polycells]
    return polycells


def _scale_cell(cell, scale=None, expand=None):
    if (scale == 1.0) or (expand == 0.0):
        return cell
    if scale is not None:
        return sscale(cell, xfact=scale, yfact=scale)
    if expand is not None:
        bbox = cell.bounds
        expanded_bbox = [
            bbox[0] - expand,
            bbox[1] - expand,
            bbox[2] + expand,
            bbox[3] + expand,
        ]
        return box(*expanded_bbox)


def _neighborshapes(polycells, scale=None, expand=None):
    nbcells = {}
    polycells = [loads(c) for c in polycells]
    for i, c in enumerate(polycells):
        c.index = i
    tree = STRtree(polycells)
    for i, cell in enumerate(polycells):
        scaled_cell = _scale_cell(cell, scale, expand)
        result = tree.query(scaled_cell)
        neighs = [n.index for n in result if n.index != i]
        nbcells[i] = neighs

    return nbcells


def _neighborpoints(cells, expand):
    cells = [eval(c) for c in cells]
    tree = cKDTree(cells)
    result = tree.query_ball_point(cells, expand)
    for i, arr in zip(range(len(cells)), result):
        try:
            arr.remove(i)
        except ValueError:
            pass
    return dict(zip(range(len(cells)), result))


if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`",
        )

    _polygonize_cells_mp = ray.remote(_polygonize_cells)
    _neighborshapes_mp = ray.remote(_neighborshapes)
    _neighborpoints_mp = ray.remote(_neighborpoints)


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
        if groupby is None:
            groupby = CONFIG.EXP_OBS
        if type_key is None:
            type_key = CONFIG.CELL_TYPE_KEY
        if centroid_key is None:
            centroid_key = CONFIG.CENTROID_KEY
        if shape_key is None:
            shape_key = CONFIG.SHAPE_KEY
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
        self.__names = [n for n, _ in self.__groups]
        self.__polycellsdb = {n: 0 for n in self.__names}
        self.__neighborsdb = {n: 0 for n in self.__names}
        self.__uniquetypes = np.unique(self.__data[type_key])
        self.__types = {n: 0 for n in self.__names}

    # @timer(prefix="Finding cell neighbors")
    def find_neighbors(
        self,
        expand: Optional[float] = None,
        scale: Optional[float] = None,
        mp: Optional[bool] = None,
    ):
        """To find the neighbors of each cell

        Args:
            expand: If the cell is shape, it means how much units to expand each cell; If the cell is point, it's the
            search radius
            scale: how much to scale each cell, only if cell is shape
            mp: whether to enable multi processing

        """
        if mp is None:
            mp = CONFIG.MULTI_PROCESSING

        # handle parameters
        if (expand is None) & (scale is None):
            raise ValueError("Neither 'expand' or 'scale' are specific")
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

        # parallel processing
        if mp & (CONFIG.OS in ["Linux", "Darwin"]):

            def exec_iterator(obj_ids):
                while obj_ids:
                    done, obj_ids = ray.wait(obj_ids)
                    yield ray.get(done[0])

            # prepare data for shape neighbor search
            if (self.__geom == "shape") & (not self.__polycells):
                results = []
                names = []
                for n, g in self.__groups:
                    results.append(_polygonize_cells_mp.remote(self.__shapekey, g))
                    names.append(n)

                    if self.__typekey is not None:
                        types = list(g[self.__typekey])
                        self.__types[n] = types

                for _ in tqdm(
                    exec_iterator(results),
                    **CONFIG.tqdm(total=len(results), desc="polygonize cells",),
                ):
                    pass

                results = ray.get(results)

                for i, n in enumerate(names):
                    self.__polycellsdb[n] = results[i]

                self.__polycells = True

            results = []
            names = []
            # shape neighbor search
            if self.__geom == "shape":
                for n, polycells in self.__polycellsdb.items():
                    results.append(_neighborshapes_mp.remote(polycells, scale, expand))
                    names.append(n)
                for _ in tqdm(
                    exec_iterator(results),
                    **CONFIG.tqdm(total=len(results), desc="find neighbors"),
                ):
                    pass
                results = ray.get(results)
            # point neighbor search
            else:
                for n, g in self.__groups:
                    results.append(_neighborpoints_mp.remote(g[self.__centkey], expand))
                    names.append(n)

                    if self.__typekey is not None:
                        types = list(g[self.__typekey])
                        self.__types[n] = types

                for _ in tqdm(
                    exec_iterator(results),
                    **CONFIG.tqdm(total=len(results), desc="find neighbors"),
                ):
                    pass
                results = ray.get(results)

            for i, n in enumerate(names):
                self.__neighborsdb[n] = results[i]
            self.__neighborsbuilt = True

        else:
            if (self.__geom == "shape") & (not self.__polycells):
                for n, g in tqdm(
                    self.__groups, **CONFIG.tqdm(desc="polygonize cells"),
                ):
                    polycells = _polygonize_cells(self.__shapekey, g)
                    self.__polycellsdb[n] = polycells

                    if self.__typekey is not None:
                        types = list(g[self.__typekey])
                        self.__types[n] = types

                self.__polycells = True
            # shape neighbor search
            if self.__geom == "shape":
                for n, polycells in tqdm(
                    self.__polycellsdb.items(), **CONFIG.tqdm(desc="find neighbors"),
                ):
                    nbcells = _neighborshapes(polycells, scale, expand)
                    self.__neighborsdb[n] = nbcells
            # point neighbor search
            else:
                for n, g in tqdm(self.__groups, **CONFIG.tqdm(desc="find neighbors"),):
                    nbcells = _neighborpoints(g[self.__centkey], expand)
                    self.__neighborsdb[n] = nbcells

                    if self.__typekey is not None:
                        types = list(g[self.__typekey])
                        self.__types[n] = types

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
        if not self.__neighborsbuilt:
            return None
        new_graphs = {n: 0 for n in self.__names}
        for n, g in self.__groups:
            centroids = [eval(c) for c in g[self.__centkey]]
            # X = [c[0] for c in centroids]
            # Y = [c[1] for c in centroids]
            vertices = [
                {"name": i, "x": x, "y": y} for i, (x, y) in enumerate(centroids)
            ]
            edges = self.__neighborsdb[n]
            graph_edges = []
            for k, vs in edges.items():
                if len(vs) > 0:
                    for v in vs:
                        distance = euclidean(centroids[k], centroids[v])
                        # graph_edges.append((k, v, distance))
                        graph_edges.append(
                            {"source": k, "target": v, "weight": distance}
                        )
                # else:
                # graph_edges.append((k, k, 0))
                # graph_edges.append({"source": k, "target": k, "weight": 0})
            # g = ig.Graph.TupleList(graph_edges, weights=True).simplify()
            g = ig.Graph.DictList(vertices, graph_edges)
            # g.vs["type"] = self.__types
            # g.vs["x"] = X
            # g.vs["y"] = Y
            new_graphs[n] = g

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
