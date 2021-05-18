from ast import literal_eval
from collections import Counter
from time import time
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .config import ANALYSIS, CONFIG
from .typing import Array
from .utils import doc
from .utils.io import df2adata_uns
from .utils.log import log_print, pretty_time


class Timer:
    """Base Class for timing"""

    _task_name: Optional[str] = None
    display_name: str
    method: Optional[str] = None
    start_time: float
    end_time: float
    used: str

    def start_timer(self) -> None:
        log_print(
            f":hourglass_not_done: [green]{self.display_name}[/green]", custom=True
        )
        if self.method is not None:
            log_print(f":hammer_and_wrench: Method: {self.method}")
        self.start_time = time()

    def stop_timer(self) -> None:
        self.end_time = time()
        self.used = pretty_time(self.end_time - self.start_time)
        log_print(f":stopwatch: [bold cyan]{self.used}[/bold cyan]", custom=True)

    @property
    def task_name(self):
        return self._task_name

    @task_name.setter
    def task_name(self, v):
        self._task_name = v
        self.display_name = ANALYSIS[v].display_name


@doc
class AnalysisBase(Timer):
    """The base class for all analysis function

    All parameters apply in this class can be used in analysis

    Args:
        data: {adata}
        task_name: The name of the analysis
        method: The method used in the run of the analysis
        exp_obs: list, How your experiments data grouped, (Default: `spatialtis.CONFIG.EXP_OBS`)
        export: bool, Whether export the result to `AnnData.uns`
        export_key: str, The name of key used to stored the exported result
        mp: bool, Whether to enable parallel processing (Default: `spatialtis.CONFIG.MULTI_PROCESSING`)

        cell_type_key: {cell_type_key}
        centroid_key: {centroid_key}
        area_key: {area_key}
        shape_key: {shape_key}
        eccentricity_key: {eccentricity_key}
        marker_key: {marker_key}
        neighbors_key: {neighbors_key}

    Attributes:
        result: To get the results
        method: The method used for the analysis, might be empty

    """

    data: AnnData
    exp_obs: List[str]
    task_name: str
    export: bool = True
    export_key: str
    mp: bool
    _result: Optional[pd.DataFrame] = None
    method: Optional[str] = None
    params: Optional[Dict] = None

    cell_type_key: str
    centroid_key: str
    area_key: str
    shape_key: str
    eccentricity_key: str
    marker_key: str
    neighbors_key: str

    def __repr__(self):
        return ""

    def __init__(
            self,
            data: AnnData,
            task_name: Optional[str] = None,
            method: Optional[str] = None,
            exp_obs: Optional[List[str]] = None,
            export: Optional[bool] = None,
            export_key: Optional[str] = None,
            mp: Optional[bool] = None,
            cell_type_key: Optional[str] = None,
            centroid_key: Optional[str] = None,
            area_key: Optional[str] = None,
            shape_key: Optional[str] = None,
            eccentricity_key: Optional[str] = None,
            marker_key: Optional[str] = None,
            neighbors_key: Optional[str] = None,
    ):
        self.data = data
        self.task_name = task_name
        if method is not None:
            self.method = method
        if exp_obs is None:
            self.exp_obs = CONFIG.EXP_OBS
            if self.exp_obs is None:
                raise ValueError("Please set CONFIG.EXP_OBS or pass `exp_obs=`")
        elif isinstance(exp_obs, (str, int, float)):
            self.exp_obs = [exp_obs]
        else:
            self.exp_obs = list(exp_obs)

        if export is not None:
            self.export = export

        if export_key is None:
            if self.task_name is not None:
                self.export_key = ANALYSIS[self.task_name].export_key
        else:
            self.export_key = export_key
        ANALYSIS[self.task_name].last_used_key = self.export_key
        if cell_type_key is None:
            self.cell_type_key = CONFIG.CELL_TYPE_KEY
        else:
            self.cell_type_key = cell_type_key
        if centroid_key is None:
            self.centroid_key = CONFIG.CENTROID_KEY
        else:
            self.centroid_key = centroid_key
        if area_key is None:
            self.area_key = CONFIG.AREA_KEY
        else:
            self.area_key = area_key
        if shape_key is None:
            self.shape_key = CONFIG.SHAPE_KEY
        else:
            self.shape_key = shape_key
        if eccentricity_key is None:
            self.eccentricity_key = CONFIG.ECCENTRICITY_KEY
        else:
            self.eccentricity_key = eccentricity_key
        if marker_key is None:
            self.marker_key = CONFIG.MARKER_KEY
        else:
            self.marker_key = marker_key
        if neighbors_key is None:
            self.neighbors_key = CONFIG.NEIGHBORS_KEY
        else:
            self.neighbors_key = neighbors_key
        if mp is None:
            self.mp = CONFIG.MP
        else:
            self.mp = mp

        if self.cell_type_key is not None:
            self.cell_types = pd.unique(self.data.obs[self.cell_type_key])

        self.neighbors_ix_key = CONFIG.neighbors_ix_key

        self.start_timer()

    def type_counter(self) -> pd.DataFrame:
        df = self.data.obs[self.exp_obs + [self.cell_type_key]]
        groups = df.groupby(self.exp_obs)
        matrix = list()
        meta = list()
        for n, g in groups:
            c = Counter(g[self.cell_type_key])
            matrix.append([c.get(t, 0) for t in self.cell_types])
            if isinstance(n, str):
                meta.append((n,))
            else:
                meta.append((*n,))
        result = dict(
            **dict(
                zip(self.exp_obs, np.asarray(meta).T),
                **dict(zip(self.cell_types, np.asarray(matrix).T)),
            )
        )
        result = pd.DataFrame(result)
        return result

    def get_exp_matrix_fraction(
            self,
            markers: Optional[Array] = None,
            types: Optional[Array] = None,
            layers_key: Optional[str] = None,
            std: Optional[float] = None,
            neighbors_ix: Optional[Array] = None,
            neighbors: Optional[tuple] = None,
            data: Optional[AnnData] = None,
    ) -> (Array, np.ndarray, AnnData):
        if data is None:
            data = self.data
        if types is not None:
            data = data[data.obs[self.cell_type_key].isin(types)].copy()
        markers_mask = []
        if markers is not None:
            if len(markers) > 1:
                markers_mask = (
                    data.var[self.marker_key].isin(markers).to_numpy(dtype=bool)
                )
            else:
                raise ValueError("Need more than two markers for `selected_markers`.")

        if std is not None:
            mask = np.asarray(data.X.std(axis=0) > std, dtype=bool)
            if len(markers_mask) == 0:
                markers_mask = mask
            else:
                markers_mask = markers_mask & mask

        if len(markers_mask) > 0:
            cut_data = data[:, markers_mask].copy()
            cut_markers = cut_data.var[self.marker_key]
        else:
            cut_data = data
            cut_markers = data.var[self.marker_key]

        if layers_key is not None:
            exp_matrix = cut_data.layers[layers_key].copy()
        else:
            exp_matrix = cut_data.X.copy()

        if neighbors is not None:
            meta = (
                cut_data.obs.reset_index(drop=True)
                    .reset_index()
                    .set_index(self.neighbors_ix_key)
            )
            cent_exp_ix = meta.loc[neighbors[0]]["index"].values
            neigh_exp_ix = meta.loc[neighbors[1]]["index"].values
            cent_exp = exp_matrix[cent_exp_ix]
            neigh_exp = exp_matrix[neigh_exp_ix]
            assert cent_exp.shape == neigh_exp.shape
            return cut_markers, (cent_exp, neigh_exp), cut_data
        elif neighbors_ix is not None:
            meta = (
                cut_data.obs.reset_index(drop=True)
                    .reset_index()
                    .set_index(self.neighbors_ix_key)
            )
            exp_ix = meta.loc[neighbors_ix]["index"].values
            exp = exp_matrix[exp_ix]
            return cut_markers, exp, cut_data
        else:
            return cut_markers, exp_matrix, cut_data

    def get_neighbors_ix(self) -> (List, List):
        need_eval = self.is_col_str(self.neighbors_key)
        if need_eval:
            neighbors = [literal_eval(n) for n in self.data.obs[self.neighbors_key]]
        else:
            neighbors = [n for n in self.data.obs[self.neighbors_key]]
        cent = [nix for nix in self.data.obs[self.neighbors_ix_key]]
        return cent, neighbors

    def get_neighbors_ix_map(self) -> Dict:
        """To get the array of index for both center and it's neighbor cells"""
        neighbors_map = {}
        cent, neighbors = self.get_neighbors_ix()
        for ix, nxs in zip(cent, neighbors):
            neighbors_map[ix] = []
            for nx in nxs:
                if ix < nx:
                    neighbors_map[ix].append(nx)

        return neighbors_map

    def get_neighbors_ix_pair(self) -> (List, List):
        """To get the array of index for both center and it's neighbor cells"""
        cent_cells = []
        neigh_cells = []
        cent, neighbors = self.get_neighbors_ix()
        for ix, nxs in zip(cent, neighbors):
            for nx in nxs:
                if ix < nx:
                    cent_cells.append(ix)
                    neigh_cells.append(nx)

        return cent_cells, neigh_cells

    def get_types_neighbors_ix(self, selected_types=None):

        if selected_types is None:
            types = self.cell_types
        else:
            types = []
            for t in selected_types:
                if t in self.cell_types:
                    types.append(t)

        neighbors = {i: {i: ([], []) for i in types} for i in types}
        # we get pairs that's not repeated from this function
        cent_cells, neigh_cells = self.get_neighbors_ix_pair()
        types_map = self.data.obs[
            [self.neighbors_ix_key, self.cell_type_key]
        ].set_index(self.neighbors_ix_key)
        cent_type = types_map.loc[cent_cells][self.cell_type_key]
        neigh_type = types_map.loc[neigh_cells][self.cell_type_key]
        for cent, neigh, c_type, n_type in zip(
                cent_cells, neigh_cells, cent_type, neigh_type
        ):
            if (c_type in types) & (n_type in types):
                # it's a pair, so we need to add it twice
                container = neighbors[c_type][n_type]
                container[0].append(cent)
                container[1].append(neigh)

                container = neighbors[n_type][c_type]
                container[0].append(neigh)
                container[1].append(cent)
        return neighbors

    def is_col_str(self, key) -> bool:
        """To determine whether a column need to eval from str

        When writing to file, the python structure like list or tuple won't be correctly interpreted,
        May need to do it manually.

        Args:
            key: The key in anndata.obs

        Returns: bool

        """
        if isinstance(self.data.obs[key][0], str):
            return True
        else:
            return False

    def export_result(self) -> None:
        export_params = {"exp_obs": self.exp_obs, "method": self.method}
        if self.params is not None:
            for k, v in self.params.items():
                export_params[k] = v
        if self.export:
            df2adata_uns(self.result, self.data, self.export_key, params=export_params)

    @property
    def neighbors_exists(self) -> bool:
        if self.neighbors_key in self.data.obs.keys():
            return True
        else:
            return False

    @property
    def result(self):
        return self._result

    @result.setter
    def result(self, v):
        self._result = v
        self.export_result()
        self.stop_timer()
