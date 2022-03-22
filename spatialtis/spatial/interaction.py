import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import CellCombs

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, options_guard
from spatialtis.utils.io import read_neighbors


@doc
def cell_interaction(data: AnnData,
                     method: str = "pval",
                     resample: int = 1000,
                     pval: float = 0.01,
                     **kwargs, ):
    """`Profiling cell-cell interaction <about/implementation.html#profiling-of-cell-cell-interaction>`_
    using permutation test

    Neighborhood analysis tells you the relationship between different type of cells

    - Association (1)
    - Avoidance (-1)
    - No relationship (0)

    Args:
        data: {adata}
        method: "pval" and "zscore" (Default: "pval")
        resample: Number of times to perform resample
        pval: {pval}
        order: If False, (Cell_A, Cell_B) and (Cell_B, Cell_A) are the same interaction (Default: False)
        **kwargs: {analysis_kwargs}

    .. seealso:: :class:`spatialtis.spatial_enrichment`

    """
    method = options_guard(method, ["pval", "zscore"])
    display_method = {"pval": "pseudo p-value", "zscore": "z-score"}
    ab = AnalysisBase(data, method=display_method[method], **kwargs)

    cc = CellCombs(ab.cell_types)

    results_data = []
    roi_tracker = []
    repeat_time = 0
    for roi_name, roi_data in ab.roi_iter(desc="Cell interaction"):
        neighbors = read_neighbors(roi_data, ab.neighbors_key)
        labels = roi_data[ab.cell_id_key]
        cell_types = roi_data[ab.cell_type_key]
        result = cc.bootstrap(
            cell_types,
            neighbors,
            labels,
            times=resample,
            pval=pval,
            method=method,
        )
        results_data += result
        roi_tracker += [roi_name]
        repeat_time = len(result)
    df = pd.DataFrame(data=results_data, columns=["type1", "type2", "value"])
    ix = pd.DataFrame(data=np.repeat(roi_tracker, repeat_time, axis=0), columns=ab.exp_obs)
    df = pd.concat([ix, df], axis=1)
    df = df.pivot_table(values="value", index=ab.exp_obs, columns=["type1", "type2"])

    ab.result = df
