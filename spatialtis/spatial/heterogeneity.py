import warnings
from ast import literal_eval
from typing import Optional, Union

import pandas as pd
from anndata import AnnData
from scipy.stats import entropy
from spatialentropy import altieri_entropy, leibovici_entropy
from tqdm import tqdm

from spatialtis.abc import AnalysisBase
from spatialtis.config import CONFIG
from spatialtis.typing import Array, Number
from spatialtis.utils import create_remote, doc, run_ray


@doc
class spatial_heterogeneity(AnalysisBase):
    """Evaluate tissue heterogeneity based on entropy

        Entropy describes the amount of information.

        - `Shannon entropy <../about/implementation.html#shannon-entropy>`_ (No spatial info included):\
            To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences\
            for each sample within the group are computed, smaller value indicates less difference within group.
        - `Leibovici entropy <../about/implementation.html#leibovici-entropy>`_:\
        You can specific the distance threshold to determine co-occurrence events.
        - `Altieri entropy <../about/implementation.html#altieri-entropy>`_:\
        You can specific the distance interval to determine co-occurrence events.

        Args:
            data: {adata}
            method: "shannon", "leibovici" and "altieri" (Default: "leibovici")
            base: The log base
            d: The distance threshold to determine co-occurrence events (method="leibovici")
            cut: Distance interval (method="altieri")
            compare: Compute Kullback-Leibler divergences based on which level (method="shannon")
            **kwargs: {analysis_kwargs}


    """

    def __init__(
        self,
        data: AnnData,
        method: str = "leibovici",
        base: Optional[Number] = None,
        d: Optional[int] = None,
        cut: Union[int, Array, None] = None,
        compare: Optional[str] = None,
        **kwargs,
    ):

        if method not in ["shannon", "altieri", "leibovici"]:
            raise ValueError(
                "Available entropy methods are 'shannon', 'altieri', 'leibovici'."
            )

        super().__init__(
            data,
            task_name="spatial_heterogeneity",
            method=f"{method.capitalize()} entropy",
            **kwargs,
        )

        if method == "shannon":
            df = self.type_counter()
            meta = df[self.exp_obs]
            df = df.iloc[:, len(self.exp_obs) : :]
            if len(df.columns) == 1:
                warnings.warn(
                    "No heterogeneity, you only have one type of cell.", UserWarning
                )
            else:
                KL_div = dict()
                ent, KL, KL_level = list(), list(), list()
                if compare is None:
                    for _, row in df.iterrows():
                        ent.append(entropy(row, base=base))
                else:
                    df = pd.concat([meta[compare], df], axis=1)
                    groups = df.groupby(compare)
                    for n, g in groups:
                        KL_div[n] = list(g.sum().values[1::])
                    for ix, row in df.iterrows():
                        compare_level = row.values[0]
                        pk = list(row.values[1::])
                        KL.append(entropy(pk, KL_div[compare_level], base=base))
                        KL_level.append(compare_level)
                        ent.append(entropy(pk, base=base))

                data = {"heterogeneity": ent}
                if compare is not None:
                    data["KL"] = KL
                    data["level"] = KL_level
                self.result = pd.concat([meta, pd.DataFrame(data)], axis=1)

        else:

            def altieri_entropy_mp(*args, **kw):
                return altieri_entropy(*args, **kw).entropy

            def leibovici_entropy_mp(*args, **kw):
                return leibovici_entropy(*args, **kw).entropy

            if method == "altieri":
                entropy_func = altieri_entropy_mp
            else:
                entropy_func = leibovici_entropy_mp

            df = data.obs[self.exp_obs + [self.cell_type_key, self.centroid_key]]

            ent = list()
            names = list()
            groups = df.groupby(self.exp_obs)
            need_eval = self.is_col_str(self.centroid_key)

            if self.mp:
                entropy_mp = create_remote(entropy_func)
                jobs = []
                for n, g in groups:
                    names.append(n)
                    types = list(g[self.cell_type_key])
                    if need_eval:
                        points = [literal_eval(i) for i in g[self.centroid_key]]
                    else:
                        points = [i for i in g[self.centroid_key]]
                    if method == "altieri":
                        jobs.append(
                            entropy_mp.remote(points, types, cut=cut, base=base)
                        )
                    else:
                        jobs.append(entropy_mp.remote(points, types, d=d, base=base))

                mp_results = run_ray(jobs, desc="Calculating heterogeneity")

                for e in mp_results:
                    ent.append(e)

            else:
                for n, g in tqdm(
                    groups, **CONFIG.pbar(desc="Calculating heterogeneity")
                ):
                    names.append(n)
                    types = list(g[self.cell_type_key])
                    if need_eval:
                        points = [literal_eval(i) for i in g[self.centroid_key]]
                    else:
                        points = [i for i in g[self.centroid_key]]
                    if method == "altieri":
                        e = altieri_entropy(points, types, cut=cut, base=base)
                    else:
                        e = leibovici_entropy(points, types, d=d, base=base)
                    ent.append(e.entropy)
            roi_heterogeneity = pd.DataFrame({"heterogeneity": ent})
            meta = pd.DataFrame(data=names, columns=self.exp_obs)
            self.result = pd.concat([meta, roi_heterogeneity], axis=1)
