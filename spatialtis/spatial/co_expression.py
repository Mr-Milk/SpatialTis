from itertools import combinations_with_replacement
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import pearsonr, spearmanr

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError, job_cutter
from spatialtis.typing import Array, Number
from spatialtis.utils import create_remote, doc, run_ray
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


def corr(
        marker_ix,
        markers,
        cent_exp,
        neigh_exp,
        corr_func,
        pval,
        corr_cutoff,
        exp_std_cutoff,
):
    result_data = []
    for p1, p2 in marker_ix:
        v1 = cent_exp[p1]
        v2 = neigh_exp[p2]
        guard1 = v1.std() > exp_std_cutoff
        guard2 = v2.std() > exp_std_cutoff
        # add a guard that make sure not ZeroDivision warning occur
        if guard1 & guard2:
            corr_value, pvalue = corr_func(v1, v2)
            if (pvalue < pval) & (abs(corr_value) > corr_cutoff):
                result_data.append([markers[p1], markers[p2], corr_value, pvalue])
    return result_data


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
            selected_types: Optional[Array] = None,
            selected_markers: Optional[Array] = None,
            layers_key: Optional[str] = None,
            exp_std_cutoff: float = 1.0,
            pval: Number = 0.01,
            corr_cutoff: Number = 0.5,
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
        super().__init__(data, task_name="spatial_co_expression", **kwargs)
        self.params = {"use_cell_type": use_cell_type}

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if self.cell_type_key not in self.data.obs.keys():
            use_cell_type = False

        if use_cell_type:
            result_data = []
            neighbors = self.get_types_neighbors_ix(selected_types)

            markers, exp_matrix, cut_data = self.get_exp_matrix_fraction(
                markers=selected_markers, layers_key=layers_key,
            )

            if self.mp:
                corr_mp = create_remote(corr_types)
                jobs = []
                for cent, neighbrs_map in neighbors.items():
                    for neigh, (cent_cells, neigh_cells) in neighbrs_map.items():
                        # to get the exp matrix of neighors and center cells
                        meta = (
                            cut_data.obs.reset_index(drop=True)
                                .reset_index()
                                .set_index(self.neighbors_ix_key)
                        )
                        cent_exp_ix = meta.loc[cent_cells]["index"].values
                        neigh_exp_ix = meta.loc[neigh_cells]["index"].values
                        cent_exp = exp_matrix[cent_exp_ix]
                        neigh_exp = exp_matrix[neigh_exp_ix]

                        jobs.append(
                            corr_mp.remote(
                                markers,
                                cent_exp.T,
                                neigh_exp.T,
                                cent,
                                neigh,
                                corr_func,
                                pval,
                                corr_cutoff,
                                exp_std_cutoff,
                            )
                        )
                results = run_ray(jobs, desc="co-expression")
                for r in results:
                    result_data += r

            else:
                for cent, neighbrs_map in pbar_iter(
                        neighbors.items(), desc="co-expression",
                ):
                    for neigh, (cent_cells, neigh_cells) in neighbrs_map.items():
                        meta = (
                            cut_data.obs.reset_index(drop=True)
                                .reset_index()
                                .set_index(self.neighbors_ix_key)
                        )
                        cent_exp_ix = meta.loc[cent_cells]["index"].values
                        neigh_exp_ix = meta.loc[neigh_cells]["index"].values
                        cent_exp = exp_matrix[cent_exp_ix]
                        neigh_exp = exp_matrix[neigh_exp_ix]
                        result_data += corr_types(
                            markers,
                            cent_exp.T,
                            neigh_exp.T,
                            cent,
                            neigh,
                            corr_func,
                            pval,
                            corr_cutoff,
                            exp_std_cutoff,
                        )

            self.result = pd.DataFrame(
                result_data,
                columns=["cell1", "marker1", "cell2", "marker2", "corr", "pvalue"],
            )

        else:
            result_data = []

            markers, exp_matrix, cut_data = self.get_exp_matrix_fraction(
                markers=selected_markers, layers_key=layers_key
            )
            meta = (
                cut_data.obs.reset_index(drop=True)
                    .reset_index()
                    .set_index(self.neighbors_ix_key)
            )
            cent_cells, neigh_cells = self.get_neighbors_ix_pair()
            cent_exp_ix = meta.loc[cent_cells]["index"].values
            neigh_exp_ix = meta.loc[neigh_cells]["index"].values

            cent_exp = exp_matrix[cent_exp_ix].T
            neigh_exp = exp_matrix[neigh_exp_ix].T
            marker_ix = np.asarray(
                [i for i in combinations_with_replacement(range(len(markers)), 2)]
            )

            if self.mp:
                corr_mp = create_remote(corr)
                jobs = []
                jobs_num = int(len(marker_ix) / 2000) + 1
                # the correlation is small computation step
                # we need to add more into one step and then parallel
                # to compensate the overhead
                for d1, d2 in job_cutter(len(marker_ix), jobs_num):
                    jobs.append(
                        corr_mp.remote(
                            marker_ix[d1:d2],
                            markers,
                            cent_exp,
                            neigh_exp,
                            corr_func,
                            pval,
                            corr_cutoff,
                            exp_std_cutoff,
                        )
                    )
                results = run_ray(jobs, desc="co-expression")
                for r in results:
                    result_data += r
            else:
                for p1, p2 in pbar_iter(marker_ix, desc="co-expression"):
                    v1 = cent_exp[p1]
                    v2 = neigh_exp[p2]
                    if (v1.sum() != 0) & (v2.sum() != 0):
                        corr_value, pvalue = corr_func(v1, v2)
                        if (pvalue < pval) & (abs(corr_value) > corr_cutoff):
                            result_data.append(
                                [markers[p1], markers[p2], corr_value, pvalue]
                            )
            self.result = pd.DataFrame(
                result_data, columns=["marker1", "marker2", "corr", "pvalue"]
            )
