import shutil
from ast import literal_eval
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from colorama import Fore

from spatialtis.config import CONFIG

from .log import logger
from .params import get_default_params


def writer_verbose(key, verbose: Optional[bool] = None):
    if verbose is None:
        verbose = CONFIG.VERBOSE.ANNDATA

    if verbose:
        if CONFIG.WORKING_ENV is not None:
            logger.info(
                f"""{Fore.GREEN}Added to AnnData, uns: {Fore.CYAN}'{key}'{Fore.RESET}"""
            )
        else:
            logger.info(f"""Added to AnnData, uns: '{key}'""")


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
    writer_verbose(key, verbose)


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
    writer_verbose(key, verbose)


def adata_uns2df(
    adata: AnnData, key: str, params: bool = False,
):
    """Read `pandas.DataFrame` object from `AnnData.uns` written by df2adata_uns

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


def filter_adata(
    adata, groupby, type_col, *keys, selected_types=None, reset_index=True
):
    """(Private) Filter `AnnData.obs` (`pandas.DataFrame`)"""

    keys = [k for k in keys if k is not None]
    df = adata.obs[groupby + keys + [type_col]]
    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]
    if reset_index:
        df = df.reset_index()

    return df


@get_default_params
def prepare_svca(
    adata: AnnData,
    export: Union[Path, str],
    entry_folder: str = "svca_data",
    groupby: Union[Sequence, str, None] = None,
    marker_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
):
    """Prepare data for SVCA analysis

    Spatial Variance Components Analysis: `SVCA <https://github.com/damienArnol/svca>`_

    The input format is separated folder for each ROI with `expressions.txt` and `positions.txt`.

    Args:
        adata: The `AnnData` object to process
        export: The directory to store the data
        entry_folder: The name of new folder to store the data
        groupby: How your experiments data grouped, (Default: `spatialtis.CONFIG.EXP_OBS`)
        marker_key: The key to store markers in `AnnData.var` (Default: `spatialtis.CONFIG.MARKER_KEY`)
                    used in the header of the expressions.txt
        centroid_key: The key to store cell centroid in `AnnData.obs` (Default: `spatialtis.CONFIG.CENTROID_KEY`)

    """

    if marker_key not in adata.var.keys():
        raise ValueError(f"{marker_key} not in anndata object's var field")
    else:
        marker_info = adata.var[marker_key]

    keys = list(groupby) + [centroid_key]

    groups = adata.obs[keys].groupby(groupby)

    p = Path(export) / entry_folder

    if p.exists():
        shutil.rmtree(p)

    p.mkdir(exist_ok=True)

    for n, g in groups:
        folder_name = "_".join([str(i) for i in n])
        folder_path = p / folder_name
        folder_path.mkdir()

        expression_path = folder_path / "expressions.txt"
        position_path = folder_path / "positions.txt"

        # export to expressions.txt
        pd.DataFrame(adata.X[[int(i) for i in g.index]], columns=marker_info).to_csv(
            expression_path, sep="\t", index=False
        )

        # export to positions.txt
        cents = []
        for c in g[centroid_key]:
            c = literal_eval(c)
            cents.append([c[0], c[1]])
        pd.DataFrame(cents).to_csv(position_path, sep="\t", index=False, header=False)
