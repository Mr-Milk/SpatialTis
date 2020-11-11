import pytest
from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.MARKER_KEY = "Markers"
CONFIG.WORKING_ENV = None


@pytest.mark.xfail
def test_colorcycle_failed():
    sp.colorcycle("V")


def test_get_colors():
    sp.get_colors(100, "Spectral")
    sp.get_colors(100, ["Spectral", "Set3"])


def test_get_linear_colors():
    sp.get_linear_colors("Spectral")
    sp.get_linear_colors(["Spectral", "Set3"])


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_cell_map():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_map(data, ROI, return_plot=True)
    sp.cell_map(data, ROI, shape_key="not_exist_key")
    sp.cell_map(data, ROI, selected_types=['12', '2'], size=(600, 600))
    sp.cell_map(data, ROI, geom="point", title="sth", palette="Spectral")
    sp.cell_map(data, ROI, geom="point", selected_types=['12', '2'], save="test.svg")


@pytest.mark.xfail
def test_cell_map_failed_param_centroied_key():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_map(data, ROI, geom="point", centroid_key="not_exist_key")


@pytest.mark.xfail
def test_cell_map_failed_param_geom():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_map(data, ROI, geom="ttt")


def test_expression_map():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.expression_map(data, ROI, selected_types=['11', '6', '8'], order=['CD20'], display=True, return_plot=True)
    sp.expression_map(data, ROI, use="scatter", expression_max=10, expression_min=0)
    sp.expression_map(data, ROI, use="scatter", expression_min=0)
    sp.expression_map(data, ROI, use="scatter", expression_max=10)


@pytest.mark.xfail
def test_expression_map_failed_param_use():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.expression_map(data, ROI, use="ttt")
