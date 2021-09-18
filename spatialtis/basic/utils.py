from collections import Counter
from typing import List

import numpy as np
import pandas as pd
from natsort import natsorted


def bbox_size(bbox):
    x = bbox[2] - bbox[0]
    y = bbox[3] - bbox[1]
    return x * y


def bbox_eccentricity(bbox) -> float:
    x = (bbox[2] - bbox[0]) / 2.0
    y = (bbox[3] - bbox[1]) / 2.0
    if x < y:
        x, y = y, x
    return np.sqrt(1.0 - y ** 2 / x ** 2)


def type_counter(data: pd.DataFrame, exp_obs: List[str], count_on: str) -> pd.DataFrame:
    """
    The return dataframe has index of exp_obs and counted types for each columns

    Args:
        data:
        exp_obs:
        count_on:

    Returns:

    """
    types = data[count_on].unique()
    types = natsorted(types)

    matrix = []
    meta = []
    for n, g in data.groupby(exp_obs, sort=False):
        c = Counter(g[count_on])
        matrix.append([c.get(t, 0) for t in types])
        if isinstance(n, (str, int, float)):
            meta.append((n,))
        else:
            meta.append((*n,))

    index = pd.MultiIndex.from_tuples(meta)
    index.names = exp_obs
    return pd.DataFrame(data=matrix, index=index, columns=types)
