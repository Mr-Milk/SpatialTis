from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import check_neighbors
from spatialtis.utils import (
    create_remote,
    df2adata_uns,
    get_default_params,
    reuse_docstring,
    run_ray,
    timer,
)


def _max_feature(x, y, **kwargs):
    reg = RandomForestRegressor(**kwargs)
    reg.fit(x, y)
    weights = reg.feature_importances_
    max_ix = np.argmax(weights)
    max_weight = weights[max_ix]
    return [max_ix, max_weight]


@timer(prefix="Finding marker expression influenced by neighbor cells")
@get_default_params
@reuse_docstring()
def exp_neighcells(
    n: Neighbors,
    std: float = 2.0,
    marker_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
    **kwargs,
):
    """Find the neighbor cells' influence on marker expression

    Random forest regressor is built to estimate the occurrence of different types of neighbor cells on
    cell's gene expression level.
    Y is the gene expression, normalize by mean.
    X is the number of different cell types.
    For example: for cell type A with gene T 700, the neighbors of cells are type A 3, type B 2, type C 5
    then Y = 700, X = (3, 2, 5)

    A reasonable std should be set, the marker expression need to have certain degree of variance.

    Args:
        n: {n}
        std: Standard deviation, threshold to filter out markers that are not variant enough
        marker_key: {marker_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}
        **kwargs: Pass to sklearn.ensemble.RandomForestRegressor (Default: random_state=0)

    """
    if export_key is None:
        export_key = CONFIG.exp_neighcell_key
    else:
        CONFIG.exp_neighcell_key = export_key

    check_neighbors(n)

    # handle kwargs pass in randomforest
    rf_kwargs = dict(random_state=0)
    for k, v in kwargs:
        rf_kwargs[k] = v

    Y = dict(zip(n.unitypes, [[] for _ in range(len(n.unitypes))]))  # gene expression
    X = dict(zip(n.unitypes, [[] for _ in range(len(n.unitypes))]))  # cell components

    adata = n.adata
    neighbors_data = n.neighbors

    for name, roi in adata.obs.groupby(n.expobs):
        neighbors = neighbors_data[name]
        type_map = dict(zip(range(len(roi)), roi[n.type_key]))

        for (center, neighs), exp in zip(neighbors.items(), adata[roi.index].X):
            t = type_map[center]
            X[t].append(Counter([type_map[i] for i in neighs if i != center]))
            Y[t].append(list(exp))

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

            for g, m in zip(y, markers):
                if np.std(g) >= std:
                    jobs.append(_max_feature_mp.remote(x, g, **rf_kwargs))
                    m_.append(m)
                    c1_.append(c1)

        mp_results = run_ray(
            jobs,
            tqdm_config=CONFIG.tqdm(
                total=len(jobs), desc="Fitting model", unit="regressor"
            ),
        )

        for (m, c1, (max_ix, max_weights)) in zip(m_, c1_, mp_results):
            if max_weights > 0:
                c2 = n.unitypes[max_ix]
                results.append([c1, m, c2, max_weights])

    else:
        with tqdm(
            **CONFIG.tqdm(
                total=len(n.unitypes) * len(markers),
                desc="Fitting model",
                unit="regressor",
            )
        ) as pbar:
            for c1 in n.unitypes:
                y = np.asarray(Y[c1]).T
                x = np.asarray(X[c1])

                for g, m in zip(y, markers):
                    pbar.update(1)
                    if np.std(g) >= std:
                        [max_ix, max_weights] = _max_feature(x, g, **rf_kwargs)
                        if max_weights > 0:
                            c2 = n.unitypes[max_ix]
                            results.append([c1, m, c2, max_weights])
            pbar.close()

    df = (
        pd.DataFrame(
            results, columns=["Cell_1", "Marker", "(Affected_by)Cell_2", "Score"]
        )
        .sort_values("Score", ascending=False)
        .reset_index(drop=True)
    )

    if export:
        df2adata_uns(df, adata, export_key, params={"std": std})

    if return_df:
        return df
