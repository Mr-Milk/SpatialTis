import os

import pytest
from anndata import read_h5ad

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.WORKING_ENV = None


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_cell_components():
    data = pytest.data
    pytest.data = data
    st.cell_components(data)


def test_cell_components_plot():
    data = pytest.data
    sp.cell_components(data, ["Patient", "Part"], save="test.svg", title="cell density", percentage=True)
    os.remove('test.svg')
    sp.cell_components(data, ["Patient", "Part"], direction="horizontal", save="test.png")
    os.remove('test.png')


def test_cell_density():
    data = pytest.data
    st.cell_density(data, (100, 100))


def test_cell_density_plot():
    data = pytest.data
    sp.cell_density(data, ["Patient"], title="cell density")
    sp.cell_density(data, ["Patient"], direction="horizontal", title="cell density")


def test_cell_co_occurrence():
    data = pytest.data
    st.cell_co_occurrence(data)


def test_cell_co_occurrence_plot():
    data = pytest.data
    sp.cell_co_occurrence(data, use="dot", display=False)
    sp.cell_co_occurrence(data, ["Patient", "Part"], use="heatmap", display=False)


def test_cell_morphology():
    data = pytest.data
    st.cell_morphology(data)


def test_cell_morphology_plot():
    data = pytest.data
    sp.cell_morphology(data, display=False)
