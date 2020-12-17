from typing import Optional

import numpy as np
import pandas as pd
from rich.progress import track

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


@timer(task_name="Finding neighbor marker dependent markers")
@get_default_params
@reuse_docstring()
def NMD_markers(
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
    """Identify neighbor markers dependent marker

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
        export_key = CONFIG.nmd_markers_key
    else:
        CONFIG.nmd_markers_key = export_key

    check_neighbors(n)

    reg_kwargs = dict(random_state=0)
    for k, v in kwargs:
        reg_kwargs[k] = v

    adata = n.adata
    neighbors_data = n.neighbors

    markers = adata.var[marker_key]
    markers_mapper = dict(zip(markers, range(len(markers))))

    if layer_key is not None:
        expmat = adata.layers[layer_key]
    else:
        expmat = adata.X

    Y = dict(zip(markers, expmat.T))  # center marker
    X = []

    for name, roi in n.groups:
        neighbors = neighbors_data[name]
        for center, neighs in enumerate(neighbors):
            ix = roi.iloc[[i for i in neighs if i != center], :].index
            if layer_key is not None:
                expmat = adata[ix].layers[layer_key]
            else:
                expmat = adata[ix].X
            X.append(np.sum(expmat, axis=0))
    X = np.asarray(X)

    if mp:
        _max_feature_mp = create_remote(_max_feature)

        jobs = []
        used_markers = []
        for m, y in Y.items():
            if np.std(y) >= expression_std_cutoff:
                jobs.append(_max_feature_mp.remote(X, y.T, method, **reg_kwargs))
                used_markers.append(m)

        mp_results = run_ray(
            jobs,
            dict(
                total=len(jobs),
                description="[green]Fitting model",
                disable=(not CONFIG.VERBOSE),
            ),
        )

        results = []
        for m, (max_ix, max_weights) in zip(used_markers, mp_results):
            if max_weights > 0:
                results.append((m, markers[max_ix], max_weights))

    else:
        results = []
        for m, y in track(
            Y.items(),
            total=len(markers),
            description="[green]Fitting model",
            disable=(not CONFIG.VERBOSE),
        ):
            if np.std(y) >= expression_std_cutoff:
                [max_ix, max_weights] = _max_feature(X, y, method, **reg_kwargs)
                if max_weights > 0:
                    results.append((m, markers[max_ix], max_weights))

    df = pd.DataFrame(
        data=results, columns=["Marker1", "Marker2", "Score",]
    ).sort_values("Score", ascending=False)

    if export:
        df2adata_uns(df, adata, export_key)

    if return_df:
        return df
