from ast import literal_eval
from collections import Counter
from time import time
from typing import Dict, List, Optional, Any

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from rich.progress import track

from spatialtis.config import Config, console
from spatialtis.typing import Array
from spatialtis.utils import doc, df2adata_uns, read_exp, log_print, pretty_time


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
    def task_name(self, v: str):
        self._task_name = v
        self.display_name = " ".join(v.split("_")).capitalize()


def neighbors_pairs(labels: List[int], neighbors: List[List[int]], duplicates: bool = False):
    p1, p2 = [], []
    if duplicates:
        for l, ns in zip(labels, neighbors):
            for n in ns:
                p1.append(l)
                p2.append(n)
    else:
        for l, ns in zip(labels, neighbors):
            for n in ns:
                if n > l:
                    p1.append(l)
                    p2.append(n)
    return np.array([p1, p2], dtype=np.int32)


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
    export_key: str
    mp: bool
    _result: Optional[pd.DataFrame] = None
    method: Optional[str] = None
    params: Optional[Dict] = None

    roi_key: str
    cell_type_key: str
    centroid_key: str
    marker_key: str
    neighbors_key: str = "cell_neighbors"
    cell_id_key: str = "cell_id"
    area_key: str = "area"
    eccentricity_key: str = "eccentricity"

    def __repr__(self):
        return ""

    def __init__(
            self,
            data: AnnData,
            method: Optional[str] = None,
            exp_obs: Optional[List[str]] = None,
            roi_key: Optional[str] = None,
            export_key: Optional[str] = None,
            cell_type_key: Optional[str] = None,
            centroid_key: Optional[str] = None,
            shape_key: Optional[str] = None,
            marker_key: Optional[str] = None,
            mp: Optional[bool] = None,
    ):
        self.data = data
        self.task_name = self.__class__.__name__
        self.method = method
        self.cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
        self.centroid_key = Config.centroid_key if centroid_key is None else centroid_key
        self.marker_key = Config.marker_key if marker_key is None else marker_key
        self.shape_key = Config.shape_key if shape_key is None else shape_key
        self.mp = Config.mp if mp is None else mp
        self.cell_types = natsorted(pd.unique(self.data.obs[self.cell_type_key]))
        self.markers = natsorted(pd.unique(self.data.var[self.marker_key]))

        if exp_obs is None:
            self.exp_obs = Config.exp_obs
            if self.exp_obs is None:
                raise ValueError("Please set Config.exp_obs or pass `exp_obs=`")
        elif isinstance(exp_obs, (str, int, float)):
            self.exp_obs = [exp_obs]
        else:
            self.exp_obs = list(exp_obs)

        if roi_key is None:
            self.roi_key = self.exp_obs[-1]
        else:
            if roi_key not in self.exp_obs:
                raise ValueError("The `roi_key` is not in your `exp_obs`")
            else:
                if self.exp_obs[-1] != roi_key:
                    exp_obs = self.exp_obs
                    exp_obs.remove(roi_key)
                    exp_obs.append(roi_key)
                    self.exp_obs = exp_obs
                self.roi_key = roi_key

        if export_key is None:
            self.export_key = self.task_name
        else:
            self.export_key = export_key

        self.start_timer()

    def roi_iter(self,
                 sort: bool = False,
                 desc: Optional[str] = None,
                 disable_pbar: bool = False,
                 ):
        if disable_pbar:
            disable = True
        else:
            disable = Config.verbose
        for roi_name, roi_data in track(self.data.obs.groupby(self.exp_obs, sort=sort),
                                        description=desc,
                                        disable=disable,
                                        console=console):
            yield roi_name, roi_data

    def roi_exp_iter(self,
                     selected_markers: Optional[List[Any]] = None,
                     layer_key: Optional[str] = None,
                     sort: bool = False,
                     desc: Optional[str] = None,
                     disable_pbar: bool = False,
                     ) -> (List, pd.DataFrame, List, np.ndarray):
        if disable_pbar:
            disable = True
        else:
            disable = Config.verbose

        markers_mask = self.data.var[self.marker_key].isin(selected_markers)
        markers = self.data.var[self.marker_key][markers_mask]
        for roi_name, roi_data in track(self.data.obs.groupby(self.exp_obs, sort=sort),
                                        description=desc,
                                        disable=disable,
                                        console=console):
            exp = read_exp(self.data[roi_data.index, markers_mask])
            yield roi_name, roi_data, markers, exp

    def type_counter(self) -> pd.DataFrame:
        matrix = []
        meta = []
        for roi_name, roi_data in self.roi_iter([self.cell_type_key], disable_pbar=True):
            c = Counter(roi_data[self.cell_type_key])
            matrix.append([c.get(t, 0) for t in self.cell_types])
            if isinstance(roi_name, (str, int, float)):
                meta.append((roi_name,))
            else:
                meta.append((*roi_name,))
        index = pd.MultiIndex.from_tuples(meta)
        index.names = self.exp_obs
        return pd.DataFrame(data=matrix, index=index, columns=self.cell_types)

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

    def export_result(self):
        export_params = {"exp_obs": self.exp_obs, "method": self.method}
        if self.params is not None:
            export_params = {**export_params, **self.params}
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
