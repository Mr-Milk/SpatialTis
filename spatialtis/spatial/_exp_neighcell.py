from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns, timer

from ._neighbors import Neighbors
from ._util import check_neighbors


def _max_feature(x, y, **kwargs):
    reg = RandomForestRegressor(**kwargs)
    reg.fit(x, y)
    weights = reg.feature_importances_
    max_ix = np.argmax(weights)
    max_weight = weights[max_ix]
    return [max_ix, max_weight]


@timer(prefix="Finding marker expression influenced by neighbor cells")
def exp_neighcells(
    n: Neighbors,
    std: float = 2.0,
    marker_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: bool = False,
    **kwargs,
):
    """Find the neighbor cell influence on marker expression

    Random forest regressor is built to estimate the occurrence of different types of neighbor cells on
    cell's gene expression level.
    Y is the gene expression, normalize by mean.
    X is the proportion of different cell types, sums up to 1.
    For example: for cell type A with gene T 700, the neighbors of cells are type A 3, type B 2, type C 5
    then Y = 700, X = (0.3, 0.2, 0.5)

    A reasonable std should be set, the marker expression need to have certain degree of variance.

    Args:
        n: a spatialtis.Neighbors object, neighbors are already computed
        std: threshold for standard deviation of marker expression
        marker_key: the key of marker in anndata.var (Default: spatialtis.CONFIG.MARKER_KEY)
        export: whether to export the result to anndata.uns
        export_key: the key used to export
        return_df: whether to return the result
        mp: whether to enable multiprocessing
        **kwargs: pass to sklearn.ensemble.RandomForestRegressor

    Returns:
        pandas.DataFrame

    """
    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY

    if export_key is None:
        export_key = CONFIG.exp_neighcell_key
    else:
        CONFIG.exp_neighcell_key = export_key

    check_neighbors(n)

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
        X[t] = (df / df.sum()).fillna(0).to_numpy()

    markers = adata.var[marker_key]

    results = (
        list()
    )  # [cell A, cell A's gene A, cell B that exert influences on cell A's gene A, weights]

    if mp:

        try:
            import ray
        except ImportError:
            raise ImportError(
                "You don't have ray installed or your OS don't support ray.",
                "Try `pip install ray` or use `mp=False`",
            )

        _max_feature_mp = ray.remote(_max_feature)

        def exec_iterator(obj_ids):
            while obj_ids:
                done, obj_ids = ray.wait(obj_ids)
                yield ray.get(done[0])

        m_ = list()
        c1_ = list()
        for c1 in n.unitypes:
            y = np.asarray(Y[c1]).T
            x = np.asarray(X[c1])

            for g, m in zip(y, markers):
                if np.std(g) >= std:
                    results.append(_max_feature_mp.remote(x, g, **kwargs))
                    m_.append(m)
                    c1_.append(c1)

        for _ in tqdm(
            exec_iterator(results), **CONFIG.tqdm(total=len(results), desc="fit model")
        ):
            pass

        mp_results = ray.get(results)

        results = list()

        for (m, c1, (max_ix, max_weights)) in zip(m_, c1_, mp_results):
            if max_weights > 0:
                c2 = n.unitypes[max_ix]
                results.append([c1, m, c2, max_weights])

    else:
        with tqdm(
            **CONFIG.tqdm(
                total=len(n.unitypes) * len(markers), desc="fit model", unit="regressor"
            )
        ) as pbar:
            for c1 in n.unitypes:
                y = np.asarray(Y[c1]).T
                x = np.asarray(X[c1])

                for g, m in zip(y, markers):
                    pbar.update(1)
                    if np.std(g) >= std:
                        [max_ix, max_weights] = _max_feature(x, g, **kwargs)
                        if max_weights > 0:
                            c2 = n.unitypes[max_ix]
                            results.append([c1, m, c2, max_weights])
            pbar.close()

    df = pd.DataFrame(results, columns=["Cell", "Marker", "Affected_by", "Score"])

    if export:
        df2adata_uns(df, adata, export_key, params={"std": std})

    if return_df:
        return df
