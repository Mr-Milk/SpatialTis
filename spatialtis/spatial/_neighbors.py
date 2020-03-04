from collections import OrderedDict
from typing import Optional, Sequence, Union

import igraph as ig

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial.distance import euclidean
from shapely.affinity import scale as sscale
from shapely.geometry import asMultiPoint, box
from shapely.strtree import STRtree

from spatialtis.config import CONFIG
from spatialtis.utils import col2adata_obs


class Neighbors(object):
    """Storage object for cell relationship

    Args:
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
        shape_col: the key name of cell shape in anndata.obs (Default: "cell_shape")
        centroid_col: anndata.obs key that store cell centroid info (Default: "centroid")

    """
    def __init__(self,
                 adata: AnnData,
                 groupby: Union[Sequence, str, None] = None,
                 type_col: Optional[str] = None,
                 shape_col: Optional[str] = 'cell_shape',
                 selected_types: Optional[Sequence] = None,
                 centroid_col: str = 'centroid'
                 ):

        # keys for query info from anndata
        if groupby is None:
            groupby = CONFIG.EXP_OBS
        if type_col is None:
            type_col = CONFIG.CELL_TYPE_COL

        keys = []
        if isinstance(groupby, str):
            keys.append(groupby)
        else:
            keys += groupby
        keys.append(shape_col)
        keys.append(centroid_col)

        if type_col is not None:
            keys.append(type_col)

        self.__adata = adata
        self.__data = adata.obs[keys]
        self.__typecol = type_col
        self.__shapecol = shape_col
        self.__centcol = centroid_col
        self.__expobs = groupby
        self.__adata_integrity = True

        if selected_types is not None:
            if type_col is None:
                raise ValueError("'type_col' not specific")
            self.__data = self.__data[self.__data[type_col].isin(selected_types)].reset_index()
            self.__adata_integrity = False

        self.__groups = self.__data.groupby(groupby)

        # define vars for storage
        self.__names = [n for n, _ in self.__groups]
        self.__polycellsdb = {n: 0 for n in self.__names}
        self.__neighborsdb = {n: 0 for n in self.__names}
        self.__neighbors_param = {'method': 0, 'units': 0}

        if type_col is not None:
            # get the cell types
            self.__uniquetypes = np.unique(self.__data[type_col])
            self.__types = {n: 0 for n in self.__names}

        # verify if neighbors and polycells built
        self._polycells = False
        self.__neighborsbuilt = False

    def _polygonize_cells(self, name, group):
        shapes = group[self.__shapecol]
        polycells = []
        for i, cell in enumerate(shapes):
            c = asMultiPoint(eval(cell))
            c.index = i
            polycells.append(c)
        self.__polycellsdb[name] = polycells

        if self.__typecol is not None:
            types = list(group[self.__typecol])
            self.__types[name] = types

    def _neighborcells(self, name, polycells, scale, expand):
        nbcells = {}
        tree = STRtree(polycells)
        for i, cell in enumerate(polycells):
            if (scale != 1.0) & (expand == 0):
                scaled_cell = sscale(cell, xfact=scale, yfact=scale)
            elif (expand != 0) & (scale == 1.0):
                bbox = cell.bounds
                expanded_bbox = [bbox[0] - expand, bbox[1] - expand, bbox[2] + expand, bbox[3] + expand]
                scaled_cell = box(*expanded_bbox)
            else:
                scaled_cell = cell
            result = tree.query(scaled_cell)
            neighs = [n.index for n in result if n.index != i]
            if len(neighs) > 0:
                nbcells[i] = neighs

        self.__neighborsdb[name] = nbcells

    def find_neighbors(self,
                       scale: float = 1,
                       expand: float = 0,
                       ):
        """To find the neighbors of each cell

        Args:
            scale: how much to scale each cell, (Default: 1, not scale)
            expand: how much units to expand each cell, (Default: 0, not expand)

        """
        # define how to enlarge cells based on user input
        use_scale = False
        use_expand = False
        unchanged = False
        if (scale < 1) | (expand < 0):
            return 'Shrink cell size is not allowed.'
        elif scale > 1:
            use_scale = True
            expand = 0
        elif expand > 0:
            use_expand = True
        else:
            unchanged = True

        # run search or not
        # if the param is the same, prevent repetitive work
        run_neighbors_search = False
        pre_method = self.__neighbors_param['method']
        pre_unit = self.__neighbors_param['units']
        if pre_method == 0:
            run_neighbors_search = True
        else:
            if (pre_method == 'scale') & (pre_unit != scale):
                run_neighbors_search = True
            elif (pre_method == 'expand') & (pre_unit != expand):
                run_neighbors_search = True

        if use_scale:
            self.__neighbors_param = {'method': 'scale', 'units': scale}
        if use_expand:
            self.__neighbors_param = {'method': 'expand', 'units': expand}
        if unchanged:
            self.__neighbors_param = {'method': 'unchanged', 'units': 0}

        if not self._polycells:
            for n, g in self.__groups:
                self._polygonize_cells(n, g)
            self._polycells = True

        if run_neighbors_search:
            for n, polycells in self.__polycellsdb.items():
                self._neighborcells(n, polycells, scale, expand)
            self.__neighborsbuilt = True

    def export_neighbors(self, export_key: str = 'cell_neighbors'):
        """Export computed neighbors

        The neighbors relationship are stored in dict, number is index of cell in each ROI
        Correspond to the order in anndata object

        This will export to anndata object's uns field

        Args:
            export_key: the key name to export

        """
        if not self.__neighborsbuilt:
            return "Please run .find_neighbors() before further analysis."

        if export_key in self.__adata.obs.keys():
            raise KeyError(f'export key "{export_key}" exists')

        self.__adata.uns[export_key] = {'data': str(self.__neighborsdb), 'param': self.__neighbors_param}

    def neighbors_count(self, export_key: str = 'neighbors_count', overwrite: bool = False):
        """Get how many neighbors for each cell

        This will write to anndata object's obs field,
        If 'selected_types' are used, this function won't work

        Args:
            export_key: the key name to export
            overwrite: if to overwrite existed key

        """

        if not self.__neighborsbuilt:
            return "Please run .find_neighbors() before further analysis."

        if self.__adata_integrity:
            counts = []
            for n, g in self.__groups:
                neighdict = self.__neighborsdb[n]
                for i in range(0, len(g)):
                    try:
                        count = len(neighdict[i])
                    except KeyError:
                        count = 0
                    counts.append(count)
            col2adata_obs(counts, self.__adata, export_key, overwrite)
        else:
            return 'Cannot write to incomplete anndata because "selected_types" are used.'

    def read_neighbors(self, read_key: str = 'cell_neighbors'):
        """Read computed neighbors from anndata

            Args:
                read_key: the key name to read

        """
        if read_key not in self.__adata.uns.keys():
            raise KeyError(f'{read_key} not exists.')
        self.__neighborsdb = eval(self.__adata.uns[read_key]['data'])
        self.__neighbors_param = self.__adata.uns[read_key]['param']

    def to_graphs(self):
        if not self.__neighborsbuilt:
            return None
        new_graphs = {n: 0 for n in self.__names}
        for n, g in self.__groups:
            centroids = [eval(c) for c in g[self.__centcol]]
            edges = self.__neighborsdb[n]
            graph_edges = []
            for k, vs in edges.items():
                for v in vs:
                    distance = euclidean(centroids[k], centroids[v])
                    graph_edges.append((str(k), str(v), distance))
            g = ig.Graph.TupleList(graph_edges, weights=True)
            g.vs["type"] = self.__types
            new_graphs[n] = g.simplify()

        return new_graphs

    @property
    def neighborsbuilt(self):
        return self.__neighborsbuilt

    @property
    def unitypes(self):
        if self.__typecol is not None:
            return self.__uniquetypes
        else:
            return None

    @property
    def types(self):
        if self.__typecol is not None:
            return self.__types
        else:
            return None

    @property
    def data(self):
        return self.__data

    @property
    def adata(self):
        return self.__adata

    @property
    def neighbors(self):
        """Return the computed neighbors relationship"""
        return self.__neighborsdb

    @property
    def expobs(self):
        return self.__expobs
