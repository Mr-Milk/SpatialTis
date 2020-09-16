import functools
import inspect
import logging
import shutil
from pathlib import Path
from time import time
from typing import Any, Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from colorama import Fore, init

from spatialtis.config import CONFIG

if CONFIG.OS == "Windows":
    init(convert=True)

logger = logging.getLogger(__name__)


def lprint(text, color="green", verbose=None):
    cmap = {"green": Fore.GREEN, "red": Fore.RED}

    if verbose is None:
        verbose = CONFIG.VERBOSE.INFO

    if CONFIG.WORKING_ENV is None:
        msg = f"{text}"
    else:
        msg = f"{cmap[color]}{text}{Fore.RESET}"
    if verbose:
        logger.info(msg)


def pretty_time(t):
    if t < 1:
        return f"{int(t * 1000)}ms"
    elif 1 <= t < 60:
        ms = t * 1000 % 1000
        return f"{int(t)}s{int(ms)}ms"
    elif 60 <= t < 3600:
        minute = t // 60
        second = t % 60
        return f"{int(minute)}m{int(second)}s"
    elif t >= 3600:
        hour = t // 3600
        minute = (t - (hour * 3600)) // 60
        second = t % 60
        return f"{int(hour)}h{int(minute)}m{int(second)}s"


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
            # to handle default parameters
            sig = inspect.signature(func)
            params_names = sig.parameters.keys()
            kw_names = kw.keys()

            if "groupby" in kw_names:
                kw["groupby"] = list(kw["groupby"])

            def _handle_kwargs(
                kwargs_writer, kwargs, default, error=True, error_message=None
            ):
                if kwargs in params_names:
                    if kwargs not in kw_names:
                        if default is not None:
                            kwargs_writer[kwargs] = default
                        else:
                            if error:
                                raise ValueError(error_message)

            # groupby
            _handle_kwargs(
                kw,
                "groupby",
                CONFIG.EXP_OBS,
                True,
                "Experiment observation unclear, set `spatialtis.CONFIG.EXP_OBS` or use argument `groupby=`",
            )
            _handle_kwargs(kw, "mp", CONFIG.MULTI_PROCESSING, False)

            # handle default key
            for a, b, c, d in zip(
                ["type_key", "centroid_key", "metric_key", "shape_key", "marker_key"],
                [
                    CONFIG.CELL_TYPE_KEY,
                    CONFIG.CENTROID_KEY,
                    CONFIG.ECCENTRICITY_KEY,
                    CONFIG.SHAPE_KEY,
                    CONFIG.MARKER_KEY,
                ],
                ["cell type", "centroid", "metric", "shape", "marker"],
                [
                    "CONFIG.CELL_TYPE_KEY",
                    "CONFIG.CENTROID_KEY",
                    "CONFIG.ECCENTRICITY_KEY",
                    "CONFIG.SHAPE_KEY",
                    "CONFIG.MARKER_KEY",
                ],
            ):

                _handle_kwargs(
                    kw, a, b, True, f"Either specific {c} key using `{a}=` or `{d}`"
                )

            ts = time()
            result = func(*args, **kw)
            te = time()
            if (suffix is not None) & verbose:
                lprint(suffix, verbose=verbose)
            if verbose:
                if CONFIG.WORKING_ENV is not None:
                    logger.info(
                        f"{Fore.GREEN}Finished! Used {Fore.CYAN}{pretty_time(te - ts)}{Fore.RESET}"
                    )
                else:
                    logger.info(f"Finished! Used {pretty_time(te - ts)}")
            return result

        return timed

    return timeit


def df2adata_uns(
    df: pd.DataFrame,
    adata: AnnData,
    key: str,
    params: Optional[Any] = None,
    verbose: Optional[bool] = None,
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
        params: add parameters
        verbose: whether to print message

    """
    container = dict(
        df=str(df.to_dict()), iname=df.index.names, colname=list(df.columns.names),
    )

    if params is not None:
        container["params"] = params

    adata.uns[key] = container

    if verbose is None:
        verbose = CONFIG.VERBOSE.ANNDATA

    if verbose:
        if CONFIG.WORKING_ENV is not None:
            logger.info(
                f"""{Fore.GREEN}Added to AnnData, uns: {Fore.CYAN}'{key}'{Fore.RESET}"""
            )
        else:
            logger.info(f"""Added to AnnData, uns: '{key}'""")


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
        if CONFIG.WORKING_ENV is not None:
            logger.info(
                f"""{Fore.GREEN}Added to AnnData, obs: {Fore.CYAN}'{key}'{Fore.RESET}"""
            )
        else:
            logger.info(f"""Added to AnnData, obs: '{key}'""")


def adata_uns2df(
    adata: AnnData, key: str, params: bool = False,
):
    """Read pandas.DataFrame object from anndata.uns if written by df2adata_uns

    Args:
        adata: the anndata object to work with
        key: which anndata.uns key you want to read from
        params: whether to return parameters

    """

    container = adata.uns[key]
    df = pd.DataFrame(eval(container["df"]))
    df.index.set_names(container["iname"], inplace=True)
    df.columns.set_names(container["colname"], inplace=True)

    if params:
        return df, container["params"]
    else:
        return df


def filter_adata(
    adata, groupby, type_col, *keys, selected_types=None, reset_index=True
):
    """(Private) Filter anndata.obs (pandas.DataFrame)"""

    keys = [k for k in keys if k is not None]
    df = adata.obs[groupby + keys + [type_col]]
    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]
    if reset_index:
        df = df.reset_index()

    return df


def prepare_svca(
    adata: AnnData,
    export: Union[Path, str],
    entry_folder: str = "svca_data",
    groupby: Union[Sequence, str, None] = None,
    marker_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
):
    """export anndata to SVCA analysis input formats

    Spatial Variance Components Analysis: `SVCA on Github <https://github.com/damienArnol/svca>`_
    The input format is separated folder for each ROI with `expressions.txt` and `positions.txt`.

    Args:
        adata: AnnData object
        export: where to store the data
        entry_folder: the folder to store the data
        marker_key: which key in anndata var is used in the header of the expressions.txt
        groupby: how your experiments grouped, (Default: spatialtis.CONFIG.EXP_OBS)
        centroid_key: anndata.obs key that store cell centroid info (Default: spatialtis.CONFIG.CENTROID_KEY)

    """
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

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
            c = eval(c)
            cents.append([c[0], c[1]])
        pd.DataFrame(cents).to_csv(position_path, sep="\t", index=False, header=False)
