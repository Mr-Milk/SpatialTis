from typing import Iterable, List, Optional

import pandas as pd
from rich.progress import track

from spatialtis.config import Config, console


def pbar_iter(obj: Iterable, desc: Optional[str] = None, **kwargs):
    for i in track(
        obj, disable=(not Config.verbose), console=console, description=desc, **kwargs
    ):
        yield i


def roi_iter(data, exp_obs: List, keys: List):
    df: pd.DataFrame = data.obs[exp_obs + keys]
    for roi_name, roi_data in df.groupby(exp_obs, sort=False):
        yield roi_name, roi_data
