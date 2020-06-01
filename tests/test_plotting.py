import pytest
from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_cell_map(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_map(data, ROI)


def test_expression_3d(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.expression_map(data, ROI, marker_col="Markers")
