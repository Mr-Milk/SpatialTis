from itertools import combinations_with_replacement
from typing import Optional, Callable

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import pearsonr, spearmanr

from spatialtis.abc import AnalysisBase, neighbors_pairs
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Number
from spatialtis.utils import create_remote, doc, run_ray
from spatialtis.utils.io import read_neighbors, read_exp
from spatialtis.utils.log import pbar_iter


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


corr_mp: Callable = create_remote(corr)
corr_t_mp: Callable = create_remote(corr_t)


DESCRIPTION = "co-expression"


@doc
class spatial_co_expression(AnalysisBase):
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
            pval: Number = 0.01,
            **kwargs,
    ):
        if method == "spearman":
            self.method = "spearman correlation"
            corr_func = spearmanr
        elif method == "pearson":
            self.method = "pearson correlation"
            corr_func = pearsonr
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
        markers = self.data.var[self.marker_key]

        if use_cell_type:
            type_pairs = self.data.obs[self.cell_type_key][pairs.ravel()].reshape(pairs.shape)
            types = pd.DataFrame(np.hstack([pairs, type_pairs]), columns=['p1', 'p2', 'c1', 'c2'])

            results_data = []
            if self.mp:
                jobs = []
                for (t1, t2), df in types.groupby(['c1', 'c2']):
                    exp1 = read_exp(self.data[df['p1'], :])
                    exp2 = read_exp(self.data[df['p2'], :])
                    end_index = len(exp1)
                    jobs.append(corr_t_mp(exp1, exp2, t1, t2, markers, corr_func, end_index, pval))

                results = run_ray(jobs, desc=DESCRIPTION)
                for r in results:
                    results_data += r

            else:
                for (t1, t2), df in types.groupby(['c1', 'c2']):
                    exp1 = read_exp(self.data[df['p1'], :])
                    exp2 = read_exp(self.data[df['p2'], :])
                    end_index = len(exp1)

                    for i in pbar_iter(range(end_index), desc=DESCRIPTION):
                        for t in range(i + 1, end_index):
                            r, p = corr_func(exp1[i], exp2[t])
                            if p < pval:
                                results_data.append([t1, markers[i], t2, markers[t], r, p])

            self.result = pd.DataFrame(results_data, columns=["cell1", "marker1",
                                                              "cell2", "marker2",
                                                              "corr", "pval"],)
        else:
            exp1 = read_exp(self.data[pairs[0], :])
            exp2 = read_exp(self.data[pairs[1], :])
            end_index = len(exp1)

            results_data = []
            if self.mp:
                jobs = []
                for i in range(end_index):
                    jobs.append(corr_mp(exp1, exp2, markers, corr_func, i, end_index, pval))

                results = run_ray(jobs, desc=DESCRIPTION)
                for r in results:
                    results_data += r

            else:
                for i in pbar_iter(range(end_index), desc=DESCRIPTION):
                    for t in range(i + 1, end_index):
                        r, p = corr_func(exp1[i], exp2[t])
                        if p < pval:
                            results_data.append([markers[i], markers[t], r, p])

            self.result = pd.DataFrame(results_data, columns=["marker1", "marker2", "corr", "pval"])

