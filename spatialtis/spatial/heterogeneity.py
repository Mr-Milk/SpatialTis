import pandas as pd
import warnings
from anndata import AnnData
from scipy.stats import entropy
from spatialtis_core import spatial_entropy

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, options_guard


@doc
def spatial_heterogeneity(data: AnnData,
                          method: str = "leibovici",
                          d: int = None,
                          cut: int = 3,
                          export_key: str = "heterogeneity",
                          **kwargs, ):
    """Evaluate tissue heterogeneity based on entropy

        Entropy describes the amount of information.

        - `Shannon entropy <../about/implementation.html#shannon-entropy>`_ (No spatial info included):\
            To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences\
            for each sample within the group are computed, smaller value indicates less difference within group.
        - `Leibovici entropy <../about/implementation.html#leibovici-entropy>`_:\
        You can specific the distance threshold to determine co-occurrence events.
        - `Altieri entropy <../about/implementation.html#altieri-entropy>`_:\
        You can specific the distance interval to determine co-occurrence events.

        Parameters
        ----------
        data : {adata}
        method : {'shannon', 'leibovici', 'altieri}, default: 'leibovici'
        d : float
            Parameters for method='leibovici',
            The distance threshold to determine co-occurrence events.
        cut : int
            Parameters for method='altieri',
            The number of distance interval to have.
        export_key : {export_key}
        **kwargs : {analysis_kwargs}

    """

    method = options_guard(method, ["shannon", "altieri", "leibovici"])

    ab = AnalysisBase(data,
                      method=f"{method.capitalize()} entropy",
                      display_name="Spatial heterogeneity",
                      export_key=export_key,
                      **kwargs)
    ab.check_cell_type()

    if method == "shannon":
        df = ab.type_counter()
        if len(df.columns) == 1:
            warnings.warn(
                "No heterogeneity, you only have one type of cell.", UserWarning
            )
        else:
            ent = [entropy(row) for _, row in df.iterrows()]
            ab.result = pd.DataFrame({"heterogeneity": ent}, index=df.index)

    else:
        points_collections = []
        types_collections = []
        track_ix = []
        # type_mapper = {t: i for i, t in enumerate(self.cell_types)}
        for roi_name, cell_types, points in ab.iter_roi(fields=['cell_type', 'centroid']):
            points_collections.append(points)
            types_collections.append(cell_types)
            track_ix.append(roi_name)

        ent = spatial_entropy(
            points_collections,
            types_collections,
            method=method,
            d=d,
            cut=cut,
            dims=ab.dimension
        )
        ab.result = pd.DataFrame(
            {"heterogeneity": ent},
            index=pd.MultiIndex.from_tuples(track_ix, names=ab.exp_obs),
        )
