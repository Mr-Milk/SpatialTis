from ast import literal_eval
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse
from spatialtis_core import (
    dumps_wkt_points,
    dumps_wkt_polygons,
    reads_wkt_points,
    reads_wkt_polygons,
)

from spatialtis.config import Config, console


def writer_verbose(key, part: str, verbose: Optional[bool] = None):
    if verbose is None:
        verbose = Config.verbose

    if verbose:
        console.print(f":package: [green]Added to AnnData, {part}: [bold cyan]'{key}'")


def df2adata_uns(
    df: pd.DataFrame,
    adata: AnnData,
    key: str,
    params: Optional[Dict] = None,
    verbose: Optional[bool] = None,
):
    """Write pandas.DataFrame with parameters to `AnnData.uns`

    The `AnnData` haven't fully support read/write of a `pandas.Dataframe` object,
    this is a temporal solution to store it in a `Dict`

    The meaning of each key:

        - **df**: The dataframe itself
        - **iname**: The name of index/MultiIndex
        - **colname**: The name of columns/MultiIndex
        - **params**: The parameters

    Args:
        df: The `pandas.DataFrame` object you want to write to the `AnnData.uns` field
        adata: The AnnData object for storage
        key: Which key in `AnnData.obs` key you want to write
        params: Add parameters
        verbose: Control the verbosity

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


def col2adata_obs(
    col: Sequence, adata: AnnData, key: str, verbose: Optional[bool] = None
):
    """Write an array to `AnnData.obs`

    Args:
        col: An array-like object that add to `AnnData.obs`
        adata: The `AnnData` object for storage
        key: Which key in `AnnData.obs` key you want to write
        verbose: Control the verbosity

    """
    adata.obs[key] = col
    writer_verbose(key, "obs", verbose=verbose)


def get_result(
    adata: AnnData,
    key: str,
    params: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict]]:
    """Read spatialtis result from `AnnData.uns` as `pandas.DataFrame` object

    To get the params, use `params=True`

    Args:
        adata: The `AnnData` object for storage
        key: Which key in `AnnData.uns` you want to read
        params: Whether to return parameters

    """

    container = adata.uns[key]
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


# def points2wkt(points: List[tuple]) -> List[str]:
#     return [wkt.dumps(Point(*coord), trim=True) for coord in points]
#
#
# def polygons2wkt(polygons: List[List[tuple]]) -> List[str]:
#     return [MultiPoint(polygon).convex_hull.wkt for polygon in polygons]


# def load_wkt(geom: List) -> List:
#     return [list(wkt.loads(data).coords) for data in geom]


def transform_points(
    data: AnnData,
    centroid_keys: Union[str, Tuple[str, str]],
    export_key: str = "centroid",
):
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

    data.obs[export_key] = dumps_wkt_points(points)
    Config.centroid_key = export_key


def transform_shapes(data: AnnData, shape_key: str, export_key: str = "cell_shape"):
    shapes = data.obs[shape_key].tolist()
    data.obs[export_key] = dumps_wkt_polygons(shapes)


def read_points(data: pd.DataFrame, centroid_key: str) -> List[Tuple[float, float]]:
    if centroid_key is None:
        raise KeyError("centroid_key is None")
    wkt_strings = data[centroid_key].tolist()
    return reads_wkt_points(wkt_strings)
    # return [list(wkt.loads(point).coords)[0] for point in data[centroid_key]]


def read_shapes(data: pd.DataFrame, shape_key: str) -> List[List[Tuple[float, float]]]:
    wkt_strings = data[shape_key].tolist()
    return reads_wkt_polygons(wkt_strings)
    # return [list(wkt.loads(shape).coords) for shape in data[shape_key]]


def read_neighbors(data: pd.DataFrame, neighbors_key: str) -> List[List[int]]:
    return [literal_eval(n) for n in data[neighbors_key]]


def read_exp(adata: AnnData, layers_key=None, dtype=None) -> np.ndarray:
    if layers_key is None:
        exp = adata.X
    else:
        exp = adata.layers[layers_key]

    if issparse(exp):
        exp = exp.toarray()

    if dtype is not None:
        return exp.T.astype(dtype)
    else:
        return exp.T  # transpose it, so every line match to a marker
