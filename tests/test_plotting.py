from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_cell_map(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ROI = data.obs[
        (data.obs["Patient"] == "HPAP005")
        & (data.obs["Part"] == "Tail")
        & (data.obs["ROI"] == "ROI1")
        ]
    sp.cell_map(ROI)
