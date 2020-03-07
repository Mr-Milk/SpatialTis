from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import _CONFIG

_CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
_CONFIG.CELL_TYPE_COL = "leiden"

data = read_h5ad("../tmp/small.h5ad")

ROI = data.obs[
    (data.obs["Patient"] == "HPAP005")
    & (data.obs["Part"] == "Tail")
    & (data.obs["ROI"] == "ROI1")
]
sp.cell_map(ROI)
