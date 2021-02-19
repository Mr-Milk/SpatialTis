from ast import literal_eval
from typing import Dict, Optional, Sequence, Tuple, Union

import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.console import console


def writer_verbose(key, part: str = "uns", verbose: Optional[bool] = None):
    if verbose is None:
        verbose = CONFIG.VERBOSE

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
    container = dict(
        df=str(df.to_dict()), iname=df.index.names, colname=list(df.columns.names),
    )

    if params is not None:
        container["params"] = params

    adata.uns[key] = container
    writer_verbose(key, verbose=verbose)


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
    adata: AnnData, key: str, params: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict]]:
    """Read spatialtis result from `AnnData.uns` as `pandas.DataFrame` object

    To get the params, use `params=True`

    Args:
        adata: The `AnnData` object for storage
        key: Which key in `AnnData.uns` you want to read
        params: Whether to return parameters

    """

    container = adata.uns[key]
    df = pd.DataFrame(literal_eval(container["df"]))
    df.index.set_names(container["iname"], inplace=True)
    df.columns.set_names(container["colname"], inplace=True)

    if params:
        return df, container["params"]
    else:
        return df
