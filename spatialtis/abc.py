from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
import warnings
from anndata import AnnData
from ast import literal_eval
from functools import cached_property
from natsort import natsorted
from rich.progress import track
from spatialtis_core import reads_wkt_points
from time import time
from typing import Any, Dict, List, Optional, Union, Sequence

from spatialtis.config import Config, console
from spatialtis.utils import df2adata_uns, doc, log_print, pretty_time, read_exp, read_shapes, default_args


class NeighborsNotFoundError(Exception):
    pass


class CellTypeNotFoundError(Exception):
    pass


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

    Parameters
    ----------
    data : {adata}
    method : str, default: None
        The method used in the run of the analysis.
    export_key : {export_key}
    display_name : str, default: None
        The name use to display the name of analysis.
    mp : bool, default: Config.mp
        Enable parallel processing, no effect since v0.5.0.
    exp_obs : {exp_obs}
    roi_key : {roi_key}
    cell_type_key : {cell_type_key}
    centroid_key : {centroid_key}
    shape_key : {shape_key}
    marker_key : {marker_key}

    """

    data: AnnData
    exp_obs: List[str]
    task_name: str
    export_key: str
    mp: bool
    _result: Optional[pd.DataFrame] = None
    method: Optional[str] = None
    params: Optional[Dict] = None
    verbose: bool = True

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
            centroid_key: Union[str, Sequence[str], None] = None,
            shape_key: Optional[str] = None,
            marker_key: Optional[str] = None,
            mp: Optional[bool] = None,
            display_name: Optional[str] = None,
            verbose: bool = True
    ):
        self.data = data
        self.dimension = 2
        self.verbose = verbose
        self.task_name = self.__class__.__name__
        if display_name is not None:
            self.display_name = display_name

        self.method = method
        self.cell_type_key = default_args(cell_type_key, Config.cell_type_key)
        self.centroid_key = default_args(centroid_key, Config.centroid_key)
        self.marker_key = default_args(marker_key, Config.marker_key)
        self.shape_key = default_args(shape_key, Config.shape_key)
        self.mp = default_args(mp, Config.mp)

        self.has_cell_type = False
        if (self.cell_type_key is not None) & (self.cell_type_key in self.data.obs_keys()):
            self.has_cell_type = True

        if exp_obs is None:
            self.exp_obs = Config.exp_obs
            if self.exp_obs is None:
                if roi_key is None:
                    raise ValueError("Please set `Config.exp_obs`/`Config.roi_key` or pass `exp_obs=`/`roi_key=`")
                else:
                    self.exp_obs = [roi_key]
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

        # assign unique id to each cell, in case of someone cut the data afterwards
        # this ensures the analysis still work with non-integrated AnnData
        if self.cell_id_key not in data.obs_keys():
            data.obs[self.cell_id_key] = [i for i in range(len(data.obs))]

        if export_key is None:
            self.export_key = self.task_name
        else:
            self.export_key = export_key
        if verbose:
            self.start_timer()

    @cached_property
    def markers(self):
        if self.marker_key is not None:
            return natsorted(pd.unique(self.data.var[self.marker_key]))
        else:
            return natsorted(pd.unique(self.data.var.index))

    def selected_markers(self, selected_markers=None):
        if selected_markers is None:
            return self.markers
        else:
            return natsorted(pd.unique(selected_markers))

    @cached_property
    def markers_col(self):
        if self.marker_key is not None:
            return self.data.var[self.marker_key]
        else:
            return self.data.var.index

    @cached_property
    def cell_types(self):
        if self.has_cell_type:
            return natsorted(pd.unique(self.data.obs[self.cell_type_key]))
        else:
            return []

    def _get_wkt_points(self, key):
        wkt_strings = self.data.obs[key].tolist()
        try:
            points = reads_wkt_points(wkt_strings)
        except Exception:
            raise IOError("If you have two columns, try `centroid_key=('cell_x', 'cell_y'). "
                          "If you store in one column, the centroid must be in wkt format, "
                          "try `spatialtis.transform_points`")
        return points

    def get_centroids(self) -> object | List:
        ckey = self.centroid_key
        # determine the type of centroid
        # by default, read 'spatial' from .obsm
        if ckey is None:
            if 'spatial' in self.data.obsm_keys():
                return self.data.obsm['spatial'].tolist()
            if 'centroid' in self.data.obs_keys():
                return self._get_wkt_points('centroid')
            else:
                raise ValueError(
                    "Spatial information not found, please set `Config.centroid_key` or pass `centroid_key=`.")

        if isinstance(ckey, str):
            if ckey in self.data.obs_keys():
                return self._get_wkt_points(ckey)
            if ckey in self.data.obsm_keys():
                return self.data.obsm[ckey].tolist()
            else:
                raise ValueError(f"The centroid key `{ckey}` not found in either `.obsm` or `.obs`")
        else:
            check = True
            for c in ckey:
                if c not in self.data.obs_keys():
                    check = False
            if check:
                return self.data.obs[list(ckey)].to_numpy().tolist()
            else:
                raise ValueError(f"The centroid keys `{ckey}` not found in `.obs`")

    # def roi_iter(
    #         self,
    #         sort: bool = False,
    #         desc: Optional[str] = None,
    #         disable_pbar: bool = False,
    # ):
    #     """Iterate through ROI with [roi_name, roi_data]
    #
    #     Args:
    #         sort: whether to sort the ROI
    #         desc: the pbar description
    #         disable_pbar: to disable pbar
    #
    #     """
    #     disable = disable_pbar if disable_pbar else not Config.progress_bar
    #
    #     for roi_name, roi_data in track(
    #             self.data.obs.groupby(self.exp_obs, sort=sort),
    #             description=f"[green]{desc}",
    #             disable=disable,
    #             console=console,
    #     ):
    #         if len(self.exp_obs) == 1:
    #             roi_name = [roi_name]
    #         yield roi_name, roi_data

    def iter_roi_data(self,
                      fields: List[str] = None,
                      ):
        iter_data = self.data.obs.copy()
        for f in fields:
            if f == 'centroid':
                points = self.get_centroids()
                if len(points[0]) == 3:
                    self.dimension = 3
                iter_data['__spatial_centroid'] = points

            if f == 'neighbors':
                self.check_neighbors()
                iter_data['__cell_neighbors'] = self.data.obsm[self.neighbors_key]

        return iter_data

    def iter_roi(self,
                 fields: List[str] = None,
                 filter_rois: List[str] = None,
                 sort: bool = False,
                 desc: str = None,
                 disable_pbar: bool = None,
                 selected_markers: List = None,
                 layer_key: str = None,
                 dtype: Any = None,
                 ):
        """A generator to iterate ROI

        Parameters
        ----------
        fields : list of str, {'centroid', 'exp', 'neighbors', 'cell_type', 'shape', 'label', 'index'}
            What fields to retrieve when iterate ROI.
        filter_rois : list of str
            The roi to be filtered.
        sort : bool
            Whether to sort ROI.
        desc : str
            The description in the progress bar.
        disable_pbar : bool
            Whether to disable progress bar.
        selected_markers : list of str
            The list of markers to be selected.
        layer_key : str
            The layer to use for expression.
        dtype :
            The datatype.

        """
        desc = default_args(desc, self.display_name)
        disable_pbar = default_args(disable_pbar, not Config.progress_bar)
        fields = default_args(fields, [])
        iter_data = self.iter_roi_data(fields)

        if 'exp' in fields:
            selected_markers = default_args(selected_markers, self.markers)
            markers_mask = self.markers_col.isin(selected_markers)
            markers = self.markers_col[markers_mask].to_numpy()

        if filter_rois is not None:
            iter_data = iter_data[iter_data[self.roi_key].isin(filter_rois)].copy()
        for roi_name, roi_data in track(
                iter_data.groupby(self.exp_obs, sort=sort),
                description=f"[green]{desc}",
                disable=disable_pbar,
                console=console,
        ):
            # pandas will show all categories in groupby even if there is no value
            if len(roi_data) == 0:
                continue
            if len(self.exp_obs) == 1:
                roi_name = [roi_name]
            yield_fields = [roi_name]

            for f in fields:
                if f == 'centroid':
                    yield_fields.append(roi_data['__spatial_centroid'].values.tolist())
                elif f == 'exp':
                    exp = read_exp(self.data[roi_data.index, markers_mask], layer_key=layer_key, dtype=dtype)
                    yield_fields.append(markers)  # ndarray
                    yield_fields.append(exp)  # ndarray
                elif f == 'neighbors':
                    neighbors = [literal_eval(n) for n in roi_data['__cell_neighbors'].values]
                    labels = roi_data[self.cell_id_key].tolist()
                    yield_fields.append(labels)  # list
                    yield_fields.append(neighbors)  # list
                elif f == 'cell_type':
                    if self.has_cell_type:
                        yield_fields.append(roi_data[self.cell_type_key].to_numpy())
                    else:
                        yield_fields.append(None)
                elif f == 'shape':
                    if self.shape_key is None:
                        yield_fields.append(None)
                    else:
                        yield_fields.append(read_shapes(roi_data, self.shape_key))
                elif f == 'label':
                    labels = roi_data[self.cell_id_key].tolist()
                    yield_fields.append(labels)
                elif f == 'index':
                    yield_fields.append(roi_data.index)

            yield yield_fields

    def type_counter(self) -> pd.DataFrame:
        self.check_cell_type()
        matrix = []
        meta = []
        for roi_name, cell_types in self.iter_roi(fields=['cell_type'], disable_pbar=True):
            c = Counter(cell_types)
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
        return self.neighbors_key in self.data.obsm_keys()

    def check_neighbors(self):
        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Neighbors not found! Run `spatialtis.find_neighbors` first.")

    def check_cell_type(self):
        if not self.has_cell_type:
            raise CellTypeNotFoundError("Cell Type not found! Please set `cell_type_key`")

    def is_rois_name_unique(self, warn=True):
        key_len = len(self.data.obs[self.roi_key].unique())
        obs_len = len([_ for _ in self.iter_roi(disable_pbar=True)])
        compare = key_len < obs_len
        if compare & warn:
            msg = "ROI selection may be incorrect, " \
                  "ROI number determined by roi_keys " \
                  "is different from exp_obs, " \
                  "use `spatialtis.make_roi_unique()` to get unique roi"
            warnings.warn(msg)
        return not compare

    @property
    def result(self):
        """Return the result of the analysis"""
        return self._result

    @result.setter
    def result(self, v):
        self._result = v
        self.export_result()
        self.stop_timer()
