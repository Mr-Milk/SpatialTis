import platform

from .config import CONFIG
from .preprocessing import read_all_ROIs, read_ROI
from .spatial import Neighbors
from .utils import df2adata_uns, adata_uns2df, prepare_svca

if platform.system() in ["Linux", "Darwin"]:
    import logging
    import ray

    ray.init(logging_level=logging.FATAL, ignore_reinit_error=True)
