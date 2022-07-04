from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from ast import literal_eval
from scipy.sparse import issparse
from spatialtis_core import (
    dumps_points_wkt,
    dumps_polygons_wkt,
    reads_wkt_points,
    reads_wkt_polygons,
)
from typing import Dict, List, Tuple

from spatialtis.config import Config, console


def writer_verbose(key, part: str, verbose: bool = None):
    if verbose is None:
        verbose = Config.verbose

    if verbose:
        console.print(f":package: [green]Added to AnnData, {part}: [bold cyan]'{key}'")


def df2adata_uns(
        df: pd.DataFrame,
        adata: AnnData,
        key: str,
        params: Dict = None,
        verbose: bool = None,
):
    """Write pandas.DataFrame with parameters to `AnnData.uns`

    The `AnnData` haven't fully support read/write of a `pandas.Dataframe` object,
    this is a temporal solution to store it in a `Dict`

    The meaning of each key:

        - **df**: The dataframe itself
        - **iname**: The name of index/MultiIndex
        - **colname**: The name of columns/MultiIndex
        - **params**: The parameters

    Parameters
    ----------
    df : pd.DataFrame
        The `pandas.DataFrame` object you want to write to the `AnnData.uns` field
    adata : AnnData
        The `AnnData` object for storage
    key : str
        key write to `.obs`
    params : bool
        Add parameters
    verbose : bool
        Control the verbosity

    """
    # To support writing with NaN
    df = df.fillna("nan")
    container = dict(
        df=str(df.to_dict()),
        iname=str(list(df.index.names)),
        colname=str(list(df.columns.names)),
    )

    if params is not None:
        container["params"] = params

    adata.uns[key] = container
    writer_verbose(key, "uns", verbose=verbose)


def col2adata(
        col: List | np.ndarray,
        adata: AnnData,
        key: str,
        slot: str = 'obs',
        verbose: bool = None
):
    """Write an array to `AnnData.obs`

    Parameters
    ----------
    col : array-like
        An array-like object that add to `AnnData.obs`
    adata : AnnData
        The `AnnData` object for storage
    key : str
        Which key in `AnnData.obs` key you want to write
    slot : str
        Which slot to write, `obs`, `obsm`, `obsp`
    verbose : bool
        Control the verbosity

    """
    getattr(adata, slot)[key] = col
    writer_verbose(key, slot, verbose=verbose)


def get_result(
        data: AnnData,
        key: str,
        params: bool = False,
) -> pd.DataFrame | Tuple[pd.DataFrame, Dict]:
    """Read spatialtis result from `AnnData.uns` as `pandas.DataFrame` object.

    To get the params, use `params=True`.

    >>> import spatialtis as st
    >>> st.get_result(data, 'cell_components')

    Parameters
    ----------
    data : AnnData
        The `AnnData` object to retrive result from.
    key : str
        Which key in `AnnData.uns` you want to read.
    params : bool, default: False
        Whether to return parameters.

    """

    container = data.uns[key]
    try:
        df = pd.DataFrame(literal_eval(container["df"]))
        df = df.replace("nan", np.nan)
        df.index.set_names(literal_eval(container["iname"]), inplace=True)
        df.columns.set_names(literal_eval(container["colname"]), inplace=True)

        if params:
            return df, container["params"]
        else:
            return df
    except KeyError:
        raise ValueError(f"Info stored in {key} is not a spatialtis result")


def wkt_points(
        data: AnnData,
        centroid_keys: str | Tuple[str, str],
        export_key: str = "centroid",
        write_config: bool = True,
):
    """Transform normal coordination in `AnnData.obs` to wkt-format.

    >>> import spatialtis as st
    >>> st.wkt_points(data, ('x', 'y'), export_key="centroid_wkt")

    Parameters
    ----------
    data : AnnData
        The `AnnData` to work on.
    centroid_keys : str or list of str
        The key or a tuple of keys that store X, Y coordination.
    export_key : str
        The key to export.
    write_config : bool, default: True
        Whether to update centroid key to global configuration.

    """
    if isinstance(centroid_keys, str):
        points = data.obs[centroid_keys].tolist()
        if isinstance(points[0], str):
            points = [literal_eval(p) for p in points]
    elif isinstance(centroid_keys, (Tuple, List)):
        points = [
            (float(x), float(y))
            for x, y in zip(data.obs[centroid_keys[0]], data.obs[centroid_keys[1]])
        ]
    else:
        raise TypeError("centroid_keys can either be str or (str, str)")

    data.obs[export_key] = dumps_points_wkt(points)
    if write_config:
        Config.centroid_key = export_key


def wkt_shapes(
        data: AnnData,
        shape_key: str,
        export_key: str = "cell_shape",
        write_config: bool = True,
):
    """Transform normal coordination in `AnnData.obs` to wkt-format.

    >>> import spatialtis as st
    >>> st.wkt_points(data, 'shape', export_key="shape_wkt")

    Parameters
    ----------
    data : AnnData
        The `AnnData` to work on.
    shape_key : str
        The key that store shape information.
    export_key : str
        The key to export.
    write_config : bool, default: True
        Whether to update shape key to global configuration.

    """
    shapes = data.obs[shape_key].tolist()
    if isinstance(shapes[0], str):
        shapes = [literal_eval(s) for s in shapes]
    data.obs[export_key] = dumps_polygons_wkt(shapes)
    if write_config:
        Config.shape_key = export_key


def read_points(data: pd.DataFrame, centroid_key: str) -> List[List[float]]:
    if centroid_key is None:
        raise KeyError("centroid_key is None")
    wkt_strings = data[centroid_key].to_numpy().tolist()
    try:
        points = reads_wkt_points(wkt_strings)
    except Exception:
        raise IOError("The points (cell coordination) must be in wkt format, "
                      "try spatialtis.transform_points")
    return points
    # return [list(wkt.loads(point).coords)[0] for point in data[centroid_key]]


def read_shapes(data: pd.DataFrame, shape_key: str) -> List[List[List[float]]]:
    if shape_key is None:
        raise KeyError("shape_key is None")
    wkt_strings = data[shape_key].tolist()
    try:
        shapes = reads_wkt_polygons(wkt_strings)
    except Exception:
        raise IOError("The shapes (cell shapes) must be in wkt format, "
                      "try spatialtis.transform_shapes")
    return shapes


def read_neighbors(data: pd.DataFrame, neighbors_key: str) -> List[List[int]]:
    """No neighbors check before read"""
    return [literal_eval(n) for n in data[neighbors_key]]


def read_exp(adata: AnnData, layer_key=None, dtype=None) -> np.ndarray:
    if layer_key is None:
        exp = adata.X
    else:
        exp = adata.layers[layer_key]

    if issparse(exp):
        exp = exp.toarray()

    if dtype is not None:
        return exp.T.astype(dtype)
    else:
        return exp.T  # transpose it, so every line match to a marker
