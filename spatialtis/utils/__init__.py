from colorama import init

from spatialtis.config import CONFIG

from .docs import reuse_docstring
from .io import adata_uns2df, col2adata_obs, df2adata_uns, filter_adata, prepare_svca
from .log import log_print, timer
from .parallel import create_remote, run_ray
from .params import get_default_params

if CONFIG.OS == "Windows":
    init(convert=True)
