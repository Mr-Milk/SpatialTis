import pytest
from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.WORKING_ENV = None


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_cell_map():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_map(data, ROI)
    sp.cell_map(data, ROI, geom="point")


def test_expression_map():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.expression_map(data, ROI, marker_key="Markers")
    sp.expression_map(data, ROI, use="scatter", marker_key="Markers")
