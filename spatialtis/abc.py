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

        self.start_timer()

    def type_counter(self) -> pd.DataFrame:
        df = self.data.obs[self.exp_obs + [self.cell_type_key]]
        groups = df.groupby(self.exp_obs)
        types = pd.unique(df[self.cell_type_key])
        matrix = list()
        meta = list()
        for i, (n, g) in enumerate(groups):
            c = Counter(g[self.cell_type_key])
            matrix.append([c.get(t, 0) for t in types])
            meta.append((*n,))
        result = dict(
            **dict(
                zip(self.exp_obs, np.asarray(meta).T),
                **dict(zip(types, np.asarray(matrix).T)),
            )
        )
        result = pd.DataFrame(result)
        return result

    def get_exp_matrix(
        self, markers: Optional[Array] = None, layers_key: Optional[str] = None
    ) -> (Array, np.ndarray, AnnData):
        if markers is not None:
            if len(markers) > 1:
                dt = self.data.T.copy()
                cut_data = dt[dt.obs[self.marker_key].isin(markers)].copy().T
                del dt
            else:
                raise ValueError("Need more than two markers for `selected_markers`.")
        else:
            cut_data = self.data
            markers = self.data.var[self.marker_key].tolist()

        if layers_key is not None:
            exp_matrix = cut_data.layers[layers_key].copy()
        else:
            exp_matrix = cut_data.X.copy()

        return markers, exp_matrix, cut_data

    def get_neighbors_ix(self) -> (List, List):
        """To get the array of index for both center and it's neighbor cells"""
        cent_cells = []
        neigh_cells = []
        need_eval = self.is_col_str(self.neighbors_key)
        for name, g in self.data.obs.reset_index(drop=True).groupby(self.exp_obs):
            real_ix = g.index
            neighbors = []
            for ix, nxs in enumerate(g[self.neighbors_key]):
                if need_eval:
                    nxs = literal_eval(nxs)
                for nx in nxs:
                    if ix > nx:
                        neighbors.append([ix, nx])
            arr = np.array(neighbors).T
            if len(arr) > 0:
                cent_cells += [int(real_ix[i]) for i in arr[0]]
                neigh_cells += [int(real_ix[i]) for i in arr[1]]
        return cent_cells, neigh_cells

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
