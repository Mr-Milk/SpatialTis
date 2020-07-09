import os

import pytest
from anndata import read_h5ad

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import CONFIG, Neighbors

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.MARKER_KEY = "Markers"
CONFIG.WORKING_ENV = None


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_spatial_dist(shared_datadir):
    data = pytest.data
    st.spatial_distribution(data, r=50, method="vmr")
    st.spatial_distribution(data, quad=(10, 10), method="quad")
    st.spatial_distribution(data, method="nns")


def test_spatial_dist_plot():
    data = pytest.data
    sp.spatial_distribution(data, ["Patient", "Part"], save="test.png")
    os.remove('test.png')
    sp.spatial_distribution(data, ["Patient", "Part"], use="heatmap", display=False)


def test_spatial_hetero():
    data = pytest.data
    st.spatial_heterogeneity(data)
    st.spatial_heterogeneity(data, method="shannon", compare="Patient")
    st.spatial_heterogeneity(data, method="leibovici", mp=True)


def test_spatial_hetero_plot():
    data = pytest.data
    sp.spatial_heterogeneity(data, ['Patient'])


def test_hotspot():
    data = pytest.data
    st.hotspot(data, grid_size=10)


def test_spatial_mp():
    data = pytest.data
    st.spatial_heterogeneity(data, mp=True)
    st.hotspot(data, grid_size=10, mp=True)
