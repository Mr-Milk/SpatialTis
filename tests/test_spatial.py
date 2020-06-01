import pytest
from anndata import read_h5ad
import os

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import CONFIG, Neighbors

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"
CONFIG.MARKER_COL = "Markers"


def test_spatial_dist(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data
    st.spatial_distribution(data, r=50)
    st.spatial_distribution(data, quad=(10, 10), method="quad")
    st.spatial_distribution(data, method="nns")


def test_spatial_dist_plot():
    data = pytest.data
    sp.spatial_distribution(data, ["Patient", "Part"], save="test.png")
    os.remove('test.png')
    sp.spatial_distribution(data, ["Patient", "Part"], method="heatmap", display=False)


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


def test_neighbors_shape_data():
    data = pytest.data
    n = Neighbors(data)
    n.find_neighbors(expand=5)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()
    pytest.n = n


def test_neighbors_point_data():
    data = pytest.data
    n = Neighbors(data, "point")
    n.find_neighbors(expand=10)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()


def test_community():
    n = pytest.n
    st.communities(n)


def test_community_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_communities_graph(data, ROI)
    sp.cell_communities_graph(data, ROI, method="static", display=False)


def test_cell_type_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_type_graph(data, ROI)


def test_neighborhood_analysis():
    n = pytest.n
    st.neighborhood_analysis(n, resample=50)


def test_neighborhood_analysis_plot():
    data = pytest.data
    sp.neighborhood_analysis(data, method="graph", save="test.png")
    os.remove('test.png')
    sp.neighborhood_analysis(data, ['Patient', 'Part'], method="heatmap", display=False)

    st.cell_components(data)
    sp.neighborhood_analysis(data, method="dot_matrix", display=False)


def test_spatial_enrichment_analysis():
    n = pytest.n
    st.spatial_enrichment_analysis(n, resample=50)


def test_spatial_enrichment_analysis_plot():
    data = pytest.data
    sp.spatial_enrichment_analysis(data, ["Patient", "Part"], display=False)


def test_exp_neighcell():
    n = pytest.n
    st.exp_neighcells(n)


def test_exp_neighcell_plot():
    data = pytest.data
    sp.exp_neighcells(data)


def test_spatial_mp():
    data = pytest.data
    n = pytest.n
    st.neighborhood_analysis(n, resample=50, mp=True)
    st.spatial_enrichment_analysis(n, resample=50, mp=True)
    st.spatial_heterogeneity(data, mp=True)
    st.hotspot(data, grid_size=10, mp=True)
    st.exp_neighcells(n, mp=True)
