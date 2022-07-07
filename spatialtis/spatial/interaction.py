import pandas as pd
from anndata import AnnData
from spatialtis_core import CellCombs

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, options_guard


@doc
def cell_interaction(data: AnnData,
                     method: str = "pval",
                     resample: int = 1000,
                     pval: float = 0.01,
                     export_key: str = "cell_interaction",
                     **kwargs, ):
    """`Profiling cell-cell interaction <about/implementation.html#profiling-of-cell-cell-interaction>`_
    using permutation test

    Neighborhood analysis tells you the relationship between different type of cells

    - Association (1)
    - Avoidance (-1)
    - No relationship (0)

    Parameters
    ----------
    data : {adata}
    method : {'pval', 'zscore'}, default: 'pval'
    resample : float, default: 1000
        Number of times to perform resample.
    pval : {pval}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    See Also
    --------
    :class:`spatialtis.spatial_enrichment`

    """
    method = options_guard(method, ["pval", "zscore"])
    display_method = {"pval": "pseudo p-value", "zscore": "z-score"}
    ab = AnalysisBase(data, method=display_method[method],
                      display_name="Cell interaction",
                      export_key=export_key,
                      **kwargs)

    cc = CellCombs(ab.cell_types)

    results_data = []
    roi_tracker = []
    repeat_time = 0
    for roi_name, cell_types, labels, neighbors in ab.iter_roi(fields=['cell_type', 'neighbors']):
        result = cc.bootstrap(
            cell_types,
            neighbors,
            labels,
            times=resample,
            pval=pval,
            method=method,
        )
        results_data += result
        repeat_time = len(result)
        roi_tracker += [roi_name for _ in range(repeat_time)]

    df = pd.DataFrame(data=results_data, columns=["type1", "type2", "value", "relationship"])
    ix = pd.DataFrame(data=roi_tracker, columns=ab.exp_obs).reset_index()
    df = df.set_index(pd.MultiIndex.from_frame(ix))
    # df = df.pivot_table(values="value", index=ab.exp_obs, columns=["type1", "type2"])

    ab.result = df
