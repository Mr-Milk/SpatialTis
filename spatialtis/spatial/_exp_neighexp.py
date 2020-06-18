from typing import Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import adata_uns2df, df2adata_uns

from ._neighbors import Neighbors
from ._util import check_neighbors


def _max_feature(x, y, **kwargs):
    reg = RandomForestRegressor(**kwargs)
    reg.fit(x, y)
    weights = reg.feature_importances_
    max_ix = np.argmax(weights)
    max_weight = weights[max_ix]
    return [max_ix, max_weight]


if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`",
        )

    _max_feature_mp = ray.remote(_max_feature)


def exp_neighexp(
    n: Neighbors,
    marker_col: Optional[str] = None,
    importance: float = 0.5,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: bool = False,
    **kwargs,
):
    if marker_col is None:
        marker_col = CONFIG.MARKER_KEY

    if export_key is None:
        export_key = CONFIG.exp_neighexp_key
    else:
        CONFIG.exp_neighexp_key = export_key

    check_neighbors(n)

    adata = n.adata
    neighbors_data = n.neighbors

    markers = adata.var[marker_col]
    markers_mapper = dict(zip(markers, range(len(markers))))

    interactions = adata_uns2df(adata, CONFIG.exp_neighcell_key)
    cc_mapper = dict()
    cexp_mapper = dict()
    X = dict(
        zip(
            [
                tuple(i[1])
                for i in interactions[["Affected_by", "Cell", "Marker"]].iterrows()
            ],
            [[] for _ in range(interactions.shape[0])],
        )
    )
    Y = dict(
        zip(
            [
                tuple(i[1])
                for i in interactions[["Affected_by", "Cell", "Marker"]].iterrows()
            ],
            [[] for _ in range(interactions.shape[0])],
        )
    )

    for cell, df in interactions.groupby("Affected_by"):
        cc_mapper[cell] = list(pd.unique(df["Cell"]))
    for cell, df in interactions.groupby("Cell"):
        cexp_mapper[cell] = list(pd.unique(df["Marker"]))

    all_cell_types = cc_mapper.keys()

    for name, roi in adata.obs.groupby(n.expobs):
        neighbors = neighbors_data[name]
        type_map = dict(zip(range(len(roi)), roi[n.type_key]))
        roi_exp = adata[roi.index].X

        for (center, neighs), exp in zip(neighbors.items(), roi_exp):

            centcell = type_map[center]
            if centcell in all_cell_types:
                selected_type = cc_mapper[centcell]
                selected_neigh = [
                    (i, type_map[i]) for i in neighs if type_map[i] in selected_type
                ]

                for (ic, c) in selected_neigh:
                    for gene in cexp_mapper[c]:
                        key = (centcell, c, gene)
                        X[key].append(list(exp))
                        Y[key].append(roi_exp[ic][markers_mapper[gene]])

    if mp & (CONFIG.OS in ["Linux", "Darwin"]):

        def exec_iterator(obj_ids):
            while obj_ids:
                done, obj_ids = ray.wait(obj_ids)
                yield ray.get(done[0])

        results = []
        combs = []
        for comb, arr in X.items():
            results.append(_max_feature_mp.remote(arr, Y[comb], **kwargs))
            combs.append(comb)

        for _ in tqdm(
            exec_iterator(results), **CONFIG.tqdm(total=len(results), desc="fit model")
        ):
            pass

        mp_results = ray.get(results)
        results = []
        for comb, (max_ix, max_weights) in zip(combs, mp_results):
            if max_weights > importance:
                results.append((markers[max_ix], *comb, max_weights))

    else:
        results = []
        for comb, arr in tqdm(X.items(), **CONFIG.tqdm(desc="fit model"),):
            [max_ix, max_weights] = _max_feature(arr, Y[comb], **kwargs)
            if max_weights > importance:
                results.append((markers[max_ix], *comb, max_weights))

    df = pd.DataFrame(
        data=results,
        columns=["Marker", "Affected_by", "Cell", "Affected_Marker", "Score"],
    )

    if export:
        df2adata_uns(df, adata, export_key)

    if return_df:
        return df
