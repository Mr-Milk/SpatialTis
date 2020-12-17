from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
from rich.progress import Progress

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import _max_feature, check_neighbors
from spatialtis.utils import (
    create_remote,
    df2adata_uns,
    get_default_params,
    reuse_docstring,
    run_ray,
    timer,
)


@timer(task_name="Finding neighbor cells dependent markers")
@get_default_params
@reuse_docstring()
def NCD_markers(
    n: Neighbors,
    expression_std_cutoff: float = 2.0,
    method: str = "xgboost",  # lasso, mutual_info
    layer_key: Optional[str] = None,
    marker_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
    **kwargs,
):
    """Identify neighbor cells dependent marker

    We used features selection methods to estimate the importance of different types of neighbor cells on
    marker expression level.
    Y is the marker expression, normalize by mean.
    X is the number of different cell types.
    For example: for cell type A with gene T 700, the neighbors of cells are type A 3, type B 2, type C 5
    then Y = 700, X = (3, 2, 5)

    A reasonable std should be set, the marker expression need to have certain degree of variance.

    Args:
        n: {n}
        expression_std_cutoff: Standard deviation, threshold to filter out markers that are not variant enough
        method: 'xgboost', 'lasso', 'mutual_info'
        layer_key: {layer_key}
        marker_key: {marker_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}
        **kwargs: Pass to different feature selection method (Default: random_state=0)

    """
    if export_key is None:
        export_key = CONFIG.ncd_markers_key
    else:
        CONFIG.ncd_markers_key = export_key

    check_neighbors(n)

    reg_kwargs = dict(random_state=0)
    for k, v in kwargs:
        reg_kwargs[k] = v

    Y = dict(zip(n.unitypes, [[] for _ in range(len(n.unitypes))]))  # gene expression
    X = dict(zip(n.unitypes, [[] for _ in range(len(n.unitypes))]))  # cell components

    adata = n.adata
    neighbors_data = n.neighbors

    for name, roi in n.groups:
        neighbors = neighbors_data[name]
        type_map = dict(zip(range(len(roi)), roi[n.type_key]))
        if layer_key is not None:
            expmat = adata[roi.index].layers[layer_key].tolist()
        else:
            expmat = adata[roi.index].X.tolist()
        for (center, neighs), exp in zip(enumerate(neighbors), expmat):
            t = type_map[center]
            X[t].append(Counter([type_map[i] for i in neighs if i != center]))
            Y[t].append(exp)

    t_cols = n.unitypes
    for t, d in X.items():
        df = pd.DataFrame(d, columns=t_cols).fillna(0)
        X[t] = df.to_numpy()

    markers = adata.var[marker_key]

    results = (
        list()
    )  # [cell A, cell A's gene A, cell B that exert influences on cell A's gene A, weights]

    if mp:

        _max_feature_mp = create_remote(_max_feature)

        jobs = []
        m_ = []
        c1_ = []
        for c1 in n.unitypes:

            y = np.asarray(Y[c1]).T
            x = np.asarray(X[c1])

            if np.std(y) >= expression_std_cutoff:
                for g, m in zip(y, markers):
                    jobs.append(_max_feature_mp.remote(x, g, method, **reg_kwargs))
                    m_.append(m)
                    c1_.append(c1)

        mp_results = run_ray(
            jobs,
            dict(
                total=len(jobs),
                description="Fitting model",
                disable=(not CONFIG.VERBOSE),
            ),
        )

        for (m, c1, (max_ix, max_weights)) in zip(m_, c1_, mp_results):
            if max_weights > 0:
                c2 = n.unitypes[max_ix]
                results.append([c1, m, c2, max_weights])

    else:
        with Progress(disable=(not CONFIG.VERBOSE)) as pbar:
            task = pbar.add_task("Fitting model", total=len(n.unitypes) * len(markers))
            for c1 in n.unitypes:
                y = np.asarray(Y[c1]).T
                x = np.asarray(X[c1])
                if np.std(y) >= expression_std_cutoff:
                    for g, m in zip(y, markers):
                        if np.std(g) >= expression_std_cutoff:
                            [max_ix, max_weights] = _max_feature(
                                x, g, method, **reg_kwargs
                            )
                            if max_weights > 0:
                                c2 = n.unitypes[max_ix]
                                results.append([c1, m, c2, max_weights])
                        pbar.update(task, advance=1)
                        pbar.refresh()
                else:
                    pbar.update(task, advance=len(markers))
                    pbar.refresh()

    df = (
        pd.DataFrame(
            results, columns=["Cell_1", "Marker", "(Affected_by)Cell_2", "Score"]
        )
        .sort_values("Score", ascending=False)
        .reset_index(drop=True)
    )

    if export:
        df2adata_uns(
            df,
            adata,
            export_key,
            params={"expression_std_cutoff": expression_std_cutoff, "method": method},
        )

    if return_df:
        return df
