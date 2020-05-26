import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# all analysis function import at entry
from .preprocessing import *
from .sta import *
from .spatial import *

from .config import CONFIG
from .utils import adata_uns2df, df2adata_uns, prepare_svca