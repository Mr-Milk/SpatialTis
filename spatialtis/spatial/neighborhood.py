from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import check_neighbors
from spatialtis.utils import df2adata_uns, get_default_params, reuse_docstring, timer


@timer(prefix="Running neighborhood analysis")
@get_default_params
@reuse_docstring()
def neighborhood_analysis(
    n: Neighbors,
    method: str = "pval",
    resample: int = 500,
    pval: float = 0.01,
    order: bool = True,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """`Profiling cell-cell interaction <about/implementation.html#profiling-of-cell-cell-interaction>`_ using permutation test

    Neighborhood analysis tells you the relationship between different type of cells

    - Association (1)
    - Avoidance (-1)
    - No relationship (0)

    This method is implemented in Rust, it executes in parallel automatically.

    Args:
        n: {n}
        method: Options are "pval" and "zscore"
        resample: Number of times to perform resample
        pval: {pval}
        order: if False, (Cell_A, Cell_B) and (Cell_B, Cell_A) are the same interaction.
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

    .. seealso:: `spatial_enrichment_analysis <#spatialtis.plotting.spatial_enrichment_analysis>`_


    """

    if export_key is None:
        export_key = CONFIG.neighborhood_analysis_key
    else:
        CONFIG.neighborhood_analysis_key = export_key

    try:
        import neighborhood_analysis as na
    except ImportError:
        raise ImportError("Package not found, try `pip install neighborhood_analysis`.")

    check_neighbors(n)
    types = n.unitypes
    cc = na.CellCombs(types, order)

    results = {}
    for name, value in tqdm(
        n.neighbors.items(), **CONFIG.tqdm(desc="Neighborhood analysis"),
    ):
        result = cc.bootstrap(
            n.types[name], value, resample, pval, method, ignore_self=True
        )
        result = {tuple(k): v for (k, v) in result}
        results[name] = result

    df = pd.DataFrame(results)
    df.index = pd.MultiIndex.from_tuples(df.index, names=("Cell type1", "Cell type2"))
    df.rename_axis(columns=n.expobs, inplace=True)

    if method == "pval":
        df = df.T.astype(int)
    else:
        df.replace(np.inf, 10.0, inplace=True)
        df.replace(-np.inf, -10.0, inplace=True)
        df = df.T

    if export:
        df2adata_uns(
            df,
            n.adata,
            export_key,
            params={"order": order, "method": method, "pval": pval},
        )

    if return_df:
        return df
