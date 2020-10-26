from typing import Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import check_neighbors
from spatialtis.utils import (
    adata_uns2df,
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


@timer(prefix="Finding marker expression influenced by neighbor markers")
@get_default_params
@reuse_docstring()
def exp_neighexp(
    n: Neighbors,
    score: float = 0.5,
    marker_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
    **kwargs,
):
    """Find the neighbor marker influence on marker expression after exp_neighcell

    Args:
        n: {n}
        score: Threshold to filter out markers are less importance
        marker_key: {marker_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}
        **kwargs: Pass to sklearn.ensemble.RandomForestRegressor (Default: random_state=0)

    """
    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY

    if export_key is None:
        export_key = CONFIG.exp_neighexp_key
    else:
        CONFIG.exp_neighexp_key = export_key

    check_neighbors(n)

    # handle kwargs pass in randomforest
    rf_kwargs = dict(random_state=0)
    for k, v in kwargs:
        rf_kwargs[k] = v

    adata = n.adata
    neighbors_data = n.neighbors

    if CONFIG.exp_neighcell_key not in adata.uns.keys():
        raise KeyError(
            f"{CONFIG.exp_neighcell_key} not found, please run spatialtis.exp_neighcells first."
        )

    markers = adata.var[marker_key]
    markers_mapper = dict(zip(markers, range(len(markers))))

    interactions = adata_uns2df(adata, CONFIG.exp_neighcell_key)
    interactions = interactions[interactions["Score"] >= score]

    cc_mapper = dict()
    cexp_mapper = dict()
    X = dict(
        zip(
            [
                tuple(i[1])
                for i in interactions[
                    ["(Affected_by)Cell_2", "Cell_1", "Marker"]
                ].iterrows()
            ],
            [[] for _ in range(interactions.shape[0])],
        )
    )
    Y = dict(
        zip(
            [
                tuple(i[1])
                for i in interactions[
                    ["(Affected_by)Cell_2", "Cell_1", "Marker"]
                ].iterrows()
            ],
            [[] for _ in range(interactions.shape[0])],
        )
    )

    for cell, df in interactions.groupby("(Affected_by)Cell_2"):
        cc_mapper[cell] = list(pd.unique(df["Cell_1"]))
    for cell, df in interactions.groupby("Cell_1"):
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
                neighs = [i for i in neighs if i != center]
                selected_neigh = [
                    (i, type_map[i]) for i in neighs if type_map[i] in selected_type
                ]

                for (ic, c) in selected_neigh:
                    for gene in cexp_mapper[c]:
                        key = (centcell, c, gene)
                        try:
                            X[key].append(list(exp))
                            Y[key].append(roi_exp[ic][markers_mapper[gene]])
                        except KeyError:
                            pass

    if mp:
        _max_feature_mp = create_remote(_max_feature)

        jobs = []
        combs = []
        for comb, arr in X.items():
            jobs.append(_max_feature_mp.remote(arr, Y[comb], **rf_kwargs))
            combs.append(comb)

        mp_results = run_ray(
            jobs,
            tqdm_config=CONFIG.tqdm(
                total=len(jobs), desc="Fitting model", unit="regressor"
            ),
        )

        results = []
        for comb, (max_ix, max_weights) in zip(combs, mp_results):
            if max_weights > 0:
                results.append((markers[max_ix], *comb, max_weights))

    else:
        results = []
        for comb, arr in tqdm(
            X.items(), **CONFIG.tqdm(desc="Fitting model", unit="regressor"),
        ):
            [max_ix, max_weights] = _max_feature(arr, Y[comb], **rf_kwargs)
            if max_weights > 0:
                results.append((markers[max_ix], *comb, max_weights))

    df = (
        pd.DataFrame(
            data=results,
            columns=[
                "Cell_2_Marker",
                "(Affected_by)Cell_2",
                "Cell_1",
                "Cell_1_Marker",
                "Score",
            ],
        )[["Cell_1", "Cell_1_Marker", "(Affected_by)Cell_2", "Cell_2_Marker", "Score"]]
        .sort_values("Score", ascending=False)
        .reset_index(drop=True)
    )

    if export:
        df2adata_uns(df, adata, export_key)

    if return_df:
        return df
