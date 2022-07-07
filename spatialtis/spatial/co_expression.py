import numpy as np
import pandas as pd
from anndata import AnnData
from ast import literal_eval
from itertools import combinations_with_replacement
from spatialtis_core import fast_corr
from typing import List, Literal

from spatialtis.abc import AnalysisBase, neighbors_pairs
from spatialtis.utils import doc, pbar_iter, options_guard
from spatialtis.utils.io import read_exp

DESCRIPTION = "co-expression"


@doc
def spatial_coexp(data: AnnData,
                  method: Literal["pearson", "spearman"] = "spearman",
                  use_cell_type: bool = False,
                  selected_markers: List[str] = None,
                  layer_key: str = None,
                  corr_thresh: float = 0.5,
                  export_key: str = "spatial_coexp",
                  **kwargs, ):
    """Identifying spatial co-expression markers using correlation

    The correlation is calculated within pairs of neighbor cells

    Parameters
    ----------
    data : {adata}
    method : {'spearman', 'pearson'}, default: 'spearman'
    use_cell_type : bool
        Whether to use cell type information.
    selected_markers : {selected_markers}
    corr_thresh : float, default: 0.5
        The minimum correlation value to store the result.
    layer_key : {layer_key}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """
    method = options_guard(method, ['spearman', 'pearson'])
    display_method = {"spearman": "spearman correlation",
                      "pearson": "pearson correlation"}
    ab = AnalysisBase(data,
                      method=display_method[method],
                      display_name="Spatial co-expression",
                      export_key=export_key,
                      **kwargs)
    ab.check_neighbors()

    if use_cell_type:
        ab.check_cell_type()

    ab.params = {"use_cell_type": use_cell_type}

    neighbors = [literal_eval(n) for n in data.obsm[ab.neighbors_key]]
    labels = data.obs[ab.cell_id_key]
    pairs = neighbors_pairs(labels, neighbors, duplicates=True)
    used_markers = ab.markers
    if selected_markers is not None:
        # sort the user input according to index in anndata to maintain order when we read exp
        order = {v: i for i, v in enumerate(ab.markers_col)}
        used_markers = sorted(selected_markers, key=lambda x: order[x[0]])
    markers_combs = [(x, y) for x, y in combinations_with_replacement(used_markers, 2)]

    if use_cell_type:
        pairs_pool = {}
        pairs_order = {}
        type_pairs = data.obs[ab.cell_type_key][pairs.ravel()] \
            .to_numpy().reshape(
            pairs.shape
        )
        for ix in range(type_pairs.shape[1]):
            i = type_pairs[:, ix]
            c = frozenset(i)
            if pairs_pool.get(c, 0) == 0:
                pairs_pool[c] = 1
                pairs_order[c] = i
            else:
                if (pairs_order[c] != i).any():
                    pairs[:, ix] = pairs[:, ix][::-1]
                    type_pairs[:, ix] = type_pairs[:, ix][::-1]
        types = pd.DataFrame(
            np.vstack([pairs, type_pairs]).T, columns=["p1", "p2", "c1", "c2"]
        )

        data_collector = []
        for (t1, t2), df in pbar_iter(types.groupby(["c1", "c2"]), desc=DESCRIPTION):
            exp1 = read_exp(data[df["p1"].to_numpy(dtype=int), :], dtype=np.float, layer_key=layer_key)
            exp2 = read_exp(data[df["p2"].to_numpy(dtype=int), :], dtype=np.float, layer_key=layer_key)

            r = fast_corr(exp1, exp2, method=method)
            d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
            d['cell1'] = t1
            d['cell2'] = t2
            d['corr'] = r
            d = d[(d['corr'] > corr_thresh) | (d['corr'] < -corr_thresh)]
            data_collector.append(d)
        d = pd.concat(data_collector)
        ab.result = d.sort_values('corr', ascending=False) \
            .reset_index(drop=True)

    else:
        exp1 = read_exp(data[pairs[0], :], dtype=np.float, layer_key=layer_key)
        exp2 = read_exp(data[pairs[1], :], dtype=np.float, layer_key=layer_key)
        r = fast_corr(exp1, exp2, method=method)
        d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
        d['corr'] = r
        if corr_thresh is not None:
            d = d[(d['corr'] > corr_thresh) | (d['corr'] < -corr_thresh)]
        ab.result = d.sort_values('corr', ascending=False).reset_index(drop=True)
