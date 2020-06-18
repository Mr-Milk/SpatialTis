import functools
import logging
import shutil
from pathlib import Path
from time import time
from typing import Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from colorama import Fore

from spatialtis.config import CONFIG

logger = logging.getLogger(__name__)


def lprint(text, color="green", verbose=None):
    cmap = {"green": Fore.GREEN, "red": Fore.RED}

    if verbose is None:
        verbose = CONFIG.VERBOSE.INFO

    if verbose:
        logger.info(f"{cmap[color]}{text}{Fore.RESET}")


def pretty_time(t):
    if t < 1:
        return f"{int(t * 1000)}ms"
    elif 1 <= t < 60:
        ms = t * 1000 % 1000
        return f"{int(t)}s{int(ms)}ms"
    elif 60 <= t < 3600:
        minute = t // 60
        second = t % 60
        return f"{minute}m{int(second)}s"
    elif t >= 3600:
        hour = t // 3600
        minute = (t - (hour * 3600)) // 60
        second = t % 60
        return f"{hour}h{minute}m{int(second)}s"


def timer(prefix=None, suffix=None, verbose=None):
    """
    Timer decorator to measure the time a function used

    Args:
        prefix: content to add at the front
        suffix: content to add at the tail
        verbose: whether to print

    Returns:

    """
    if verbose is None:
        verbose = CONFIG.VERBOSE.INFO

    def timeit(func):
        @functools.wraps(func)
        def timed(*args, **kw):
            if (prefix is not None) & verbose:
                lprint(prefix, verbose=verbose)
            ts = time()
            result = func(*args, **kw)
            te = time()
            if (suffix is not None) & verbose:
                lprint(suffix, verbose=verbose)
            if verbose:
                logger.info(
                    f"{Fore.GREEN}Finished! Used {Fore.CYAN}{pretty_time(te - ts)}{Fore.RESET}"
                )
            return result

        return timed

    return timeit


def df2adata_uns(
    df: pd.DataFrame, adata: AnnData, key: str, verbose: Optional[bool] = None
):
    """Preserve all info in pd.DataFrame as dict, and write to anndata.uns
    The anndata object haven't fully support read/write of a pandas.Dataframe object,
    this is an solution to store all the information in a dicts

    The meaning of each key:
    'df' The dataframe itself
    'iname' The name of index/MultiIndex
    'colname' The name of columns/MultiIndex

    Args:
        df: the pandas.DataFrame object you want to write to the anndata.uns field
        adata: the anndata object to work with
        key: which anndata.uns key you want to write to
        verbose: whether to print message

    """
    container = dict(
        df=str(df.to_dict()), iname=df.index.names, colname=list(df.columns.names),
    )

    adata.uns[key] = container

    if verbose is None:
        verbose = CONFIG.VERBOSE.ANNDATA

    if verbose:
        logger.info(
            f"""{Fore.GREEN}Added to AnnData, uns: {Fore.CYAN}'{key}'{Fore.RESET}"""
        )


def col2adata_obs(
    col: Sequence, adata: AnnData, key: str, verbose: Optional[bool] = None
):
    """Write a Sequence to anndata.obs

    Args:
        col: the Sequence object
        adata: the anndata object to work with
        key: which anndata.obs key you want to write to
        verbose: whether to print message

    """
    adata.obs[key] = col

    if verbose is None:
        verbose = CONFIG.VERBOSE.ANNDATA

    if verbose:
        logger.info(
            f"""{Fore.GREEN}Added to AnnData, obs: {Fore.CYAN}'{key}'{Fore.RESET}"""
        )


def adata_uns2df(
    adata: AnnData, key: str,
):
    """Read pandas.DataFrame object from anndata.uns if written by df2adata_uns

    Args:
        adata: the anndata object to work with
        key: which anndata.uns key you want to read from

    """

    container = adata.uns[key]
    df = pd.DataFrame(eval(container["df"]))
    df.index.set_names(container["iname"], inplace=True)
    df.columns.set_names(container["colname"], inplace=True)

    return df


def filter_adata(
    adata, groupby, type_col, *keys, selected_types=None, reset_index=True
):
    """(Private) Filter anndata.obs (pandas.DataFrame)
    """

    keys = [k for k in keys if k is not None]
    df = adata.obs[groupby + keys + [type_col]]
    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]
    if reset_index:
        df = df.reset_index()

    return df


def prepare_svca(
    adata: AnnData,
    export: Union[Path, str] = "/",
    marker_col: str = "Markers",
    *,
    groupby: Union[Sequence, str, None] = None,
    centroid_col: str = "centroid",
    entry_folder: str = "svca_data",
    overwrite: bool = True,
):
    """export anndata to SVCA analysis input formats

    Spatial Variance Components Analysis: `SVCA on Github <https://github.com/damienArnol/svca>`_
    The input format is separated folder for each ROI with `expressions.txt` and `positions.txt`.

    Args:
        adata: AnnData object
        export: where to store the data
        marker_col: which key in anndata var is used in the header of the expressions.txt
        groupby:
        centroid_col:
        entry_folder:
        overwrite: whether overwrite the export folder when already exists

    """
    if marker_col not in adata.var.keys():
        raise ValueError(f"{marker_col} not in anndata object's var field")
    else:
        marker_info = adata.var[marker_col]

    if groupby is None:
        groupby = CONFIG.EXP_OBS

    keys = list(groupby) + [centroid_col]

    groups = adata.obs[keys].groupby(groupby)

    p = Path(export) / entry_folder

    if p.exists():
        if overwrite:
            shutil.rmtree(p)
        else:
            raise FileExistsError(f"Folder `{entry_folder}` is already exists")

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
        for c in g[centroid_col]:
            c = eval(c)
            cents.append([c[0], c[1]])
        pd.DataFrame(cents).to_csv(position_path, sep="\t", index=False, header=False)
