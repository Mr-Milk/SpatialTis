import logging
import sys

from .config import CONFIG
from .preprocessing import *
from .spatial import (
    Neighbors,
    communities,
    exp_neighcells,
    exp_neighexp,
    hotspot,
    neighborhood_analysis,
    spatial_distribution,
    spatial_enrichment_analysis,
    spatial_heterogeneity,
)
from .stats import cell_co_occurrence, cell_components, cell_density, cell_morphology
from .utils import adata_uns2df, df2adata_uns, prepare_svca

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console = logging.StreamHandler(sys.stdout)
console.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(message)s")
console.setFormatter(formatter)
logger.addHandler(console)
