from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from milkviz import dot

from spatialtis import get_result
from spatialtis.utils import doc


@doc
def spatial_heterogeneity(
    data: AnnData,
    groupby: Optional[str] = None,
    key: str = "spatial_heterogeneity",
    **plot_options,
):
    pdata = get_result(data, key).reset_index()
    exp_obs = pdata.columns[0:-1].tolist()
    if (len(exp_obs) == 1) | (groupby == exp_obs[-1]):
        options = dict(color="#5DAC81", **plot_options)
        ax = sns.barplot(data=pdata, y="heterogeneity", x=exp_obs[-1], **options)
    else:
        groupby = exp_obs[0] if groupby is None else groupby
        ax = sns.boxplot(data=pdata, y="heterogeneity", x=groupby, **plot_options)
    return ax


@doc
def cell_dispersion(
    data: AnnData,
    key: str = "cell_dispersion",
    **plot_options,
):
    pdata, param = get_result(data, key, params=True)
    pdata = pd.pivot(
        pdata, columns="cell_type", index=param["exp_obs"], values="pattern"
    )
    pattern = {}
    for t, arr in pdata.iteritems():
        pattern[t] = {0: 0, 1: 0, 2: 0, 3: 0, **Counter(arr)}
    pdata = pd.DataFrame(pattern).T[[0, 1, 2, 3]]
    colors = np.repeat(
        [["#FFC408", "#c54a52", "#4a89b9", "#5a539d"]], len(pdata), axis=0
    )
    return dot(
        dot_size=pdata.to_numpy(dtype=int),
        dot_hue=colors,
        legend_title="ROI",
        xticklabels=["No Cell", "Random", "Regular", "Cluster"],
        yticklabels=pdata.index,
    )
