import pytest
from anndata import read_h5ad

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.MARKER_KEY = "Markers"
CONFIG.WORKING_ENV = None


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_spatial_dist():
    data = pytest.data
    st.spatial_distribution(data, r=50, method="vmr", export_key=CONFIG.spatial_distribution_key)
    st.spatial_distribution(data, quad=(10, 10), method="quad", return_df=True)
    st.spatial_distribution(data, method="nns")


@pytest.mark.xfail
def test_spatial_dist_failed_param_method():
    data = pytest.data
    st.spatial_distribution(data, method="ttt")


def test_spatial_dist_plot():
    data = pytest.data
    plot_kwargs = dict(
        annotated=True,
        title="sth",
        save="test.png",
        return_plot=True
    )
    sp.spatial_distribution(data, ["Patient", "Part"], **plot_kwargs)
    sp.spatial_distribution(data, ["Patient", "Part"], use="heatmap", display=False)


def test_spatial_hetero():
    data = pytest.data
    st.spatial_heterogeneity(data, method="shannon", compare="Patient")
    st.spatial_heterogeneity(data, method="leibovici")
    st.spatial_heterogeneity(data, method="altieri")


@pytest.mark.xfail
def test_spatial_hetero_failed_param_method():
    data = pytest.data
    st.spatial_heterogeneity(data, method="what")


def test_spatial_hetero_plot():
    data = pytest.data
    sp.spatial_heterogeneity(data, ['Patient'])


def test_hotspot():
    data = pytest.data
    st.hotspot(data, grid_size=10)


def test_spatial_dist_mp():
    data = pytest.data
    st.spatial_distribution(data, mp=True)


def test_spatial_hetero_mp():
    data = pytest.data
    st.spatial_heterogeneity(data, method="leibovici", mp=True)
    st.spatial_heterogeneity(data, method="altieri", mp=True)


def test_hotspot_mp():
    data = pytest.data
    st.hotspot(data, grid_size=10, mp=True)
