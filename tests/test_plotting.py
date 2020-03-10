import pytest
from anndata import read_h5ad
from pathlib import Path

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"

DATA_DIR = str(Path("tests/data/small.h5ad").absolute())


def test_cell_map():
    data = read_h5ad(DATA_DIR)
    ROI = data.obs[
        (data.obs["Patient"] == "HPAP005")
        & (data.obs["Part"] == "Tail")
        & (data.obs["ROI"] == "ROI1")
        ]
    sp.cell_map(ROI)
