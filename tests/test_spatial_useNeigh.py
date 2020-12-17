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
    n.find_neighbors()
    n.find_neighbors(scale=2)
    n.find_neighbors(expand=3)
    n.find_neighbors(expand=3)
    n.find_neighbors(scale=2, expand=3)
    n.neighbors_count(export_key=CONFIG.neighbors_count_key)
    CONFIG.WORKING_ENV = "jupyter_notebook"
    n.export_neighbors(export_key=CONFIG.neighbors_key)
    CONFIG.WORKING_ENV = "None"
    n.read_neighbors()
    pytest.n = n


def test_neighbors_point_data():
    data = pytest.data
    n = Neighbors(data, "point")
    n.find_neighbors(expand=6)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()


@pytest.mark.xfail
def test_neighbors_failed_param_geom():
    data = pytest.data
    n = Neighbors(data, "ttt")


@pytest.mark.xfail
def test_neighbors_failed_point_no_expand():
    data = pytest.data
    n = Neighbors(data, "point")
    n.find_neighbors()


@pytest.mark.xfail
def test_neighbors_failed_expand_lt_0():
    data = pytest.data
    n = Neighbors(data)
    n.find_neighbors(expand=-15)


@pytest.mark.xfail
def test_neighbors_failed_scale_lt_1():
    data = pytest.data
    n = Neighbors(data)
    n.find_neighbors(scale=0.5)


@pytest.mark.xfail
def test_neighbors_failed_neighbors_not_run():
    data = pytest.data
    n = Neighbors(data)
    n.export_neighbors()


@pytest.mark.xfail
def test_neighbors_failed_neighbors_not_count():
    data = pytest.data
    n = Neighbors(data)
    n.neighbors_count()


def test_neighbors_plot():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_neighbors(data, ROI, save="test.png")
    sp.cell_neighbors(data, ROI, use="_static", display=False, save="test.png")


def test_community():
    n = pytest.n
    st.communities(n)
    st.communities(n, export_key=CONFIG.community_key, return_df=True)


def test_community_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_communities(data, ROI)
    sp.cell_communities(data, ROI, use="_static", display=False)


def test_cell_type_graph():
    data = pytest.data
    ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
    sp.cell_neighbors(data, ROI)


def test_neighborhood_analysis():
    n = pytest.n
    st.neighborhood_analysis(n, resample=50)
    st.neighborhood_analysis(n, resample=50, method="zscore", export_key='zscore_na', return_df=True)
    st.neighborhood_analysis(n, resample=50, order=False, export_key="unorder_na")
    st.neighborhood_analysis(n, resample=50, order=False, method="zscore", export_key='zscore_unorder_na')


def test_neighborhood_analysis_plot():
    data = pytest.data
    sp.neighborhood_analysis(data, ['Patient', 'Part'], use="heatmap", display=False)
    sp.neighborhood_analysis(data, ['Patient', 'Part'], use="heatmap", display=False, key='zscore_na')
    sp.neighborhood_analysis(data, ['Patient', 'Part'], use="heatmap", display=False, key='unorder_na')
    sp.neighborhood_analysis(data, ['Patient', 'Part'], use="heatmap", display=False, key='zscore_unorder_na')

    st.cell_components(data)
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, save="test.png")
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, key="zscore_na")
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, key="unorder_na")
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, key="unorder_na")
    sp.neighborhood_analysis(data, use="dot_matrix", display=False, key="zscore_unorder_na")


def test_spatial_enrichment_analysis():
    n = pytest.n
    st.spatial_enrichment_analysis(n, threshold=0.1, selected_markers=['CD20', 'CD57'], resample=50)
    data = n.adata
    data.layers['new_layer'] = data.X >= 0.1
    st.spatial_enrichment_analysis(n, threshold=0.1, layers_key='new_layer', resample=50,
                                   export_key=CONFIG.spatial_enrichment_analysis_key,
                                   return_df=True)


@pytest.mark.xfail
def test_spatial_enrichment_analysis_failed_param_threshold():
    n = pytest.n
    st.spatial_enrichment_analysis(n, resample=50)


@pytest.mark.xfail
def test_spatial_enrichment_analysis_failed_param_selected_markers():
    n = pytest.n
    st.spatial_enrichment_analysis(n, threshold=0.1, selected_markers=['CD20'], resample=50)


def test_spatial_enrichment_analysis_plot():
    data = pytest.data
    sp.spatial_enrichment_analysis(data, display=False)


def test_ncd_markers():
    n = pytest.n
    st.NCD_markers(n)


def test_ncd_markers_plot():
    data = pytest.data
    sp.NCD_markers(data)


def test_nmd_markers():
    n = pytest.n
    st.NMD_markers(n)


def test_spatial_mp():
    n = pytest.n
    st.exp_neighcells(n, mp=True)
    st.nmd_markers(n, mp=True)
