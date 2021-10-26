from ast import literal_eval
from collections import Counter
from time import time
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from rich.progress import track

from spatialtis.config import Config, console
from spatialtis.typing import Array
from spatialtis.utils import df2adata_uns, doc, log_print, pretty_time, read_exp


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


def neighbors_pairs(
    labels: List[int], neighbors: List[List[int]], duplicates: bool = False
):
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
        method: The method used in the run of the analysis
        export_key: {export_key}
        display_name: The name use to display the name of analysis
        mp: bool, Enable parallel processing (Default: :code:`spatialtis.Config.mp`), not apply to most of the analysis
        exp_obs: {exp_obs}
        roi_key: {roi_key}
        cell_type_key: {cell_type_key}
        centroid_key: {centroid_key}
        shape_key: {shape_key}
        marker_key: {marker_key}

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
        display_name: Optional[str] = None,
    ):
        self.data = data
        self.task_name = self.__class__.__name__
        if display_name is not None:
            self.display_name = display_name

        self.method = method
        self.cell_type_key = (
            Config.cell_type_key if cell_type_key is None else cell_type_key
        )
        self.centroid_key = (
            Config.centroid_key if centroid_key is None else centroid_key
        )
        self.marker_key = Config.marker_key if marker_key is None else marker_key
        self.shape_key = Config.shape_key if shape_key is None else shape_key
        self.mp = Config.mp if mp is None else mp
        if self.cell_type_key is not None:
            self.cell_types = natsorted(pd.unique(self.data.obs[self.cell_type_key]))
        if self.marker_key is not None:
            self.markers = natsorted(pd.unique(self.data.var[self.marker_key]))
            self.markers_col = self.data.var[self.marker_key]
        else:
            self.markers = natsorted(pd.unique(self.data.var.index))
            self.markers_col = self.data.var.index

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

    def roi_iter(
        self,
        sort: bool = False,
        desc: Optional[str] = None,
        disable_pbar: bool = False,
    ):
        """Iterate through ROI with [roi_name, roi_data]

        Args:
            sort: whether to sort the ROI
            desc: the pbar description
            disable_pbar: to disable pbar

        """
        if disable_pbar:
            disable = True
        else:
            disable = not Config.progress_bar

        if len(self.exp_obs) == 1:
            for roi_name, roi_data in track(
                self.data.obs.groupby(self.exp_obs, sort=sort),
                description=f"[green]{desc}",
                disable=disable,
                console=console,
            ):
                yield [roi_name], roi_data
        else:
            for roi_name, roi_data in track(
                self.data.obs.groupby(self.exp_obs, sort=sort),
                description=f"[green]{desc}",
                disable=disable,
                console=console,
            ):
                yield roi_name, roi_data

    def roi_exp_iter(
        self,
        selected_markers: Optional[List[Any]] = None,
        layer_key: Optional[str] = None,
        dtype: Any = None,
        sort: bool = False,
        desc: Optional[str] = None,
        disable_pbar: bool = False,
    ) -> (List, pd.DataFrame, List, np.ndarray):
        if disable_pbar:
            disable = True
        else:
            disable = not Config.progress_bar
        selected_markers = (
            self.markers if selected_markers is None else selected_markers
        )
        markers_mask = self.markers_col.isin(selected_markers)
        markers = self.markers_col[markers_mask]
        if len(self.exp_obs) == 1:
            for roi_name, roi_data in track(
                self.data.obs.groupby(self.exp_obs, sort=sort),
                description=f"[green]{desc}",
                disable=disable,
                console=console,
            ):
                exp = read_exp(self.data[roi_data.index, markers_mask], layer_key=layer_key, dtype=dtype)
                yield [roi_name], roi_data, markers, exp
        else:
            for roi_name, roi_data in track(
                self.data.obs.groupby(self.exp_obs, sort=sort),
                description=f"[green]{desc}",
                disable=disable,
                console=console,
            ):
                exp = read_exp(self.data[roi_data.index, markers_mask], layer_key=layer_key, dtype=dtype)
                yield roi_name, roi_data, markers, exp

    def type_counter(self) -> pd.DataFrame:
        matrix = []
        meta = []
        for roi_name, roi_data in self.roi_iter(
            [self.cell_type_key], disable_pbar=True
        ):
            c = Counter(roi_data[self.cell_type_key])
            matrix.append([c.get(t, 0) for t in self.cell_types])
            if isinstance(roi_name, (str, int, float)):
                meta.append((roi_name,))
            else:
                meta.append((*roi_name,))
        index = pd.MultiIndex.from_tuples(meta)
        index.names = self.exp_obs
        return pd.DataFrame(data=matrix, index=index, columns=self.cell_types)

    def export_result(self):
        export_params = {"exp_obs": self.exp_obs, "method": self.method}
        if self.params is not None:
            export_params = {**export_params, **self.params}
        df2adata_uns(self.result, self.data, self.export_key, params=export_params)

    @property
    def neighbors_exists(self) -> bool:
        if self.neighbors_key in self.data.obs.keys():
            return True
        else:
            return False

    @property
    def result(self):
        """Return the result of the analysis"""
        return self._result

    @result.setter
    def result(self, v):
        self._result = v
        self.export_result()
        self.stop_timer()
