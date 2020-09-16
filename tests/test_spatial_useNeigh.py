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


def test_neighbors_shape_data():
    data = pytest.data
    n = Neighbors(data)
    n.find_neighbors(scale=2)
    n.find_neighbors(expand=3)
    n.find_neighbors(expand=3)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()
    pytest.n = n


def test_neighbors_point_data():
    data = pytest.data
    n = Neighbors(data, "point")
    n.find_neighbors(expand=6)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()


def test_neighbors_plot():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_neighbors(data, ROI)
    sp.cell_neighbors(data, ROI, use="static", display=False)


def test_community():
    n = pytest.n
    st.communities(n)


def test_community_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_communities(data, ROI)
    sp.cell_communities(data, ROI, use="static", display=False)


def test_cell_type_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_neighbors(data, ROI)


def test_neighborhood_analysis():
    n = pytest.n
    st.neighborhood_analysis(n, resample=50)
    st.neighborhood_analysis(n, resample=50, order=False, export_key="unorder_na")


def test_neighborhood_analysis_plot():
    data = pytest.data
    sp.neighborhood_analysis(data, ['Patient', 'Part'], use="heatmap", display=False)

    st.cell_components(data)
    sp.neighborhood_analysis(data, use="dot_matrix", display=False)
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, key="unorder_na")


def test_spatial_enrichment_analysis():
    n = pytest.n
    st.spatial_enrichment_analysis(n, threshold=0.1, resample=50)


def test_spatial_enrichment_analysis_plot():
    data = pytest.data
    sp.spatial_enrichment_analysis(data, display=False)


def test_exp_neighcell():
    n = pytest.n
    st.exp_neighcells(n)


def test_exp_neighcell_plot():
    data = pytest.data
    sp.exp_neighcells(data)


def test_exp_neighexp():
    n = pytest.n
    st.exp_neighexp(n)


def test_spatial_mp():
    n = pytest.n
    data = pytest.data
    n_mp = Neighbors(data)
    n_mp.find_neighbors(expand=3)
    st.neighborhood_analysis(n, resample=50)
    st.spatial_enrichment_analysis(n, threshold=0.1, resample=50)
    st.exp_neighcells(n, mp=True)
    st.exp_neighexp(n, mp=True)
