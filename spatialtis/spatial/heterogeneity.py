import warnings
from typing import Optional, Union

import pandas as pd
from anndata import AnnData
from scipy.stats import entropy
from spatialtis_core import spatial_entropy

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Array, Number
from spatialtis.utils import doc, read_points


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
            d: :code:`method="leibovici"`, The distance threshold to determine co-occurrence events
            cut: :code:`method="altieri"`, Distance interval
            **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        method: str = "leibovici",
        d: Optional[int] = None,
        cut: Union[int, Array, None] = None,
        order: bool = False,
        **kwargs,
    ):

        methods_list = ["shannon", "altieri", "leibovici"]
        if method not in methods_list:
            raise ValueError(
                f"Unknonw method: {method}," f"Available: {', '.join(methods_list)}."
            )

        super().__init__(
            data,
            method=f"{method.capitalize()} entropy",
            **kwargs,
        )

        if method == "shannon":
            df = self.type_counter()
            if len(df.columns) == 1:
                warnings.warn(
                    "No heterogeneity, you only have one type of cell.", UserWarning
                )
            else:
                ent = [entropy(row) for _, row in df.iterrows()]
                self.result = pd.DataFrame({"heterogeneity": ent}, index=df.index)

        else:
            points_collections = []
            types_collections = []
            track_ix = []
            type_mapper = {t: i for i, t in enumerate(self.cell_types)}
            for roi_name, roi_data in self.roi_iter(desc="Spatial heterogeneity"):
                points_collections.append(read_points(roi_data, self.centroid_key))
                types_collections.append(roi_data[self.cell_type_key].map(type_mapper))
                track_ix.append(roi_name)

            ent = spatial_entropy(
                points_collections,
                types_collections,
                method=method,
                d=d,
                cut=cut,
                order=order,
            )
            self.result = pd.DataFrame(
                {"heterogeneity": ent},
                index=pd.MultiIndex.from_tuples(track_ix, names=self.exp_obs),
            )
