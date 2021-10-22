from itertools import combinations_with_replacement
from typing import Callable, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import pearsonr, spearmanr
from spatialtis_core import fast_corr

from spatialtis.abc import AnalysisBase, neighbors_pairs
from spatialtis.typing import Number
from spatialtis.utils import create_remote, doc, pbar_iter, run_ray, NeighborsNotFoundError
from spatialtis.utils.io import read_exp, read_neighbors


def corr_types(
        markers,
        cent_exp,
        neigh_exp,
        cent,
        neigh,
        corr_func,
        pval,
        corr_cutoff,
        exp_std_cutoff,
):
    result_data = []
    markers_num = len(markers)
    marker_ix = np.asarray(
        [i for i in combinations_with_replacement(range(markers_num), 2)]
    )
    for p1, p2 in marker_ix:
        v1 = cent_exp[p1]
        v2 = neigh_exp[p2]
        guard1 = v1.std() > exp_std_cutoff
        guard2 = v2.std() > exp_std_cutoff
        # add a guard that make sure not ZeroDivision warning occur
        if guard1 & guard2:
            corr_value, pvalue = corr_func(v1, v2)
            if (pvalue < pval) & (abs(corr_value) > corr_cutoff):
                result_data.append(
                    [cent, markers[p1], neigh, markers[p2], corr_value, pvalue]
                )
    return result_data


def corr(exp1, exp2, markers, corr_func, start, end, pval):
    results = []
    for t in range(start + 1, end):
        r, p = corr_func(exp1[start], exp2[t])
        if p < pval:
            results.append([markers[start], markers[t], r, p])


def corr_t(exp1, exp2, t1, t2, markers, corr_func, end, pval):
    results = []
    for i in range(end):
        for t in range(i + 1, end):
            r, p = corr_func(exp1[i], exp2[t])
            if p < pval:
                results.append([t1, markers[i], t2, markers[t], r, p])


DESCRIPTION = "co-expression"


@doc
class spatial_coexp(AnalysisBase):
    """Identifying spatial co-expression markers using correlation

    The correlation is calculated within pairs of neighbor cells

    Args:
        data: {adata}
        method: "spearman" or "pearson" (Default: "spearman")
        selected_markers: {selected_markers}
        layers_key: {layers_key}
        pval: {pval}

    """

    def __init__(
            self,
            data: AnnData,
            method: str = "spearman",
            use_cell_type: bool = False,
            layer_key: Optional[str] = None,
            thresh: Number = 0.5,
            **kwargs,
    ):
        if method == "spearman":
            self.method = "spearman correlation"
        elif method == "pearson":
            self.method = "pearson correlation"
        else:
            raise ValueError("Available options are `spearman` and `pearson`.")
        super().__init__(data, **kwargs)
        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if self.cell_type_key not in self.data.obs.keys():
            use_cell_type = False

        self.params = {"use_cell_type": use_cell_type}

        neighbors = read_neighbors(self.data.obs, self.neighbors_key)
        labels = self.data.obs[self.cell_id_key]
        pairs = neighbors_pairs(labels, neighbors)
        markers = self.markers_col
        markers_combs = [(x, y) for x, y in combinations_with_replacement(markers, 2)]

        if use_cell_type:
            pairs_pool = {}
            pairs_order = {}
            type_pairs = self.data.obs[self.cell_type_key][pairs.ravel()]\
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

            results_data = []

            # if self.mp:
            #     # init remote function when needed
            #     corr_t_mp: Callable = create_remote(corr_t)
            #     jobs = []
            #     for (t1, t2), df in types.groupby(["c1", "c2"]):
            #         exp1 = read_exp(self.data[df["p1"].to_numpy(dtype=int), :])
            #         exp2 = read_exp(self.data[df["p2"].to_numpy(dtype=int), :])
            #         end_index = len(exp1)
            #         jobs.append(
            #             corr_t_mp.remote(
            #                 exp1, exp2, t1, t2, markers, corr_func, end_index, pval
            #             )
            #         )
            #     results = run_ray(jobs, desc=DESCRIPTION)
            #     for r in results:
            #         results_data += r

            # else:
            # for (t1, t2), df in pbar_iter(types.groupby(["c1", "c2"]), desc=DESCRIPTION):
            #     exp1 = read_exp(self.data[df["p1"].to_numpy(dtype=int), :])
            #     exp2 = read_exp(self.data[df["p2"].to_numpy(dtype=int), :])
            #     end_index = len(exp1)
            #     for i in range(end_index):
            #         for t in range(i + 1, end_index):
            #             r, p = corr_func(exp1[i], exp2[t])
            #             if p < pval:
            #                 results_data.append(
            #                     [t1, markers[i], t2, markers[t], r, p]
            #                 )
            #
            # self.result = pd.DataFrame(
            #     results_data,
            #     columns=["cell1", "marker1", "cell2", "marker2", "corr", "pval"],
            # )
            data_collector = []
            for (t1, t2), df in pbar_iter(types.groupby(["c1", "c2"]), desc=DESCRIPTION):
                exp1 = read_exp(self.data[df["p1"].to_numpy(dtype=int), :], dtype=np.float)
                exp2 = read_exp(self.data[df["p2"].to_numpy(dtype=int), :], dtype=np.float)

                r = fast_corr(exp1, exp2, method=method)
                d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
                d['cell1'] = t1
                d['cell2'] = t2
                d['corr'] = r
                d = d[(d['corr'] > thresh) | (d['corr'] < -thresh)]
                data_collector.append(d)
            d = pd.concat(data_collector)
            self.result = d.sort_values('corr', ascending=False) \
                .reset_index(drop=True)

        else:
            exp1 = read_exp(self.data[pairs[0], :], dtype=np.float)
            exp2 = read_exp(self.data[pairs[1], :], dtype=np.float)
            r = fast_corr(exp1, exp2, method=method)
            d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
            d['corr'] = r
            self.result = d[(d['corr'] > thresh) | (d['corr'] < -thresh)] \
                .sort_values('corr', ascending=False) \
                .reset_index(drop=True)
            #
            # results_data = []
            # if self.mp:
            #     # init remote function when needed
            #     corr_mp: Callable = create_remote(corr)
            #     jobs = []
            #     for i in range(end_index):
            #         jobs.append(
            #             corr_mp.remote(exp1, exp2, markers, corr_func, i, end_index, pval)
            #         )
            #     results = run_ray(jobs, desc=DESCRIPTION)
            #     for r in results:
            #         results_data += r
            #
            # else:
            #     for i in pbar_iter(range(end_index), desc=DESCRIPTION):
            #         for t in range(i + 1, end_index):
            #             r, p = corr_func(exp1[i], exp2[t])
            #             if p < pval:
            #                 results_data.append([markers[i], markers[t], r, p])
            #
            # self.result = pd.DataFrame(
            #     results_data, columns=["marker1", "marker2", "corr", "pval"]
            # )
