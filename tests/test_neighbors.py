import pytest

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

type_order = ["a", "b", "c"]


@pytest.mark.parametrize("method,kw", [
    ("kdtree", dict(r=30)),
    ("kdtree", dict(k=5)),
    ("kdtree", dict(k=5, r=30)),
    ("delaunay", dict())
])
def test_find_neighbors_2d(data2d, method, kw):
    Config.loads(data2d)
    st.find_neighbors(data2d, method=method, **kw)


@pytest.mark.parametrize("method,kw", [
    ("kdtree", dict(r=30)),
    ("kdtree", dict(k=5)),
    ("kdtree", dict(k=5, r=30)),
])
def test_find_neighbors_3d(data3d, method, kw):
    Config.loads(data3d)
    st.find_neighbors(data3d, method=method, **kw)


def test_find_neighbors_shape(data_shape):
    Config.loads(data_shape)
    st.find_neighbors(data_shape, scale=1.2, method="rtree", shape_key="cell_shape")


def test_community(data2d):
    Config.loads(data2d)
    st.cell_community(data2d)


@pytest.mark.parametrize("method", ["moran_i", "geary_c"])
def test_spatial_autocorr(data2d, data3d, method):
    Config.loads(data2d)
    st.spatial_autocorr(data2d, method=method)
    Config.loads(data3d)
    st.spatial_autocorr(data3d, method=method)


def test_neighborhood(data2d, data3d):
    Config.loads(data2d)
    st.cell_interaction(data2d)
    Config.loads(data3d)
    st.cell_interaction(data3d, method="zscore")


def test_neighborhood_plot(data2d):
    Config.loads(data2d)
    sp.cell_interaction(data2d)
    sp.cell_interaction(data2d, use="heatmap", type_order=type_order)
    sp.cell_interaction(data2d, use="heatmap", plot_value="value")


def test_spatial_enrichment(data2d, data3d):
    Config.loads(data2d)
    st.spatial_enrichment(data2d)
    Config.loads(data3d)
    st.spatial_enrichment(data3d, threshold=0.1)


def test_spatial_enrichment_plot(data2d):
    Config.loads(data2d)
    sp.spatial_enrichment(data2d)


def test_spatial_coexp(data2d, data3d):
    Config.loads(data2d)
    st.spatial_coexp(data2d)
    st.spatial_coexp(data2d, use_cell_type=True)

    Config.loads(data3d)
    st.spatial_coexp(data3d)
    st.spatial_coexp(data3d, use_cell_type=True)


def test_ncd_markers(data2d, data3d):
    Config.loads(data2d)
    st.NCD_marker(data2d)

    Config.loads(data3d)
    st.NCD_marker(data3d)


def test_nmd_markers(data2d, data3d):
    Config.loads(data2d)
    st.NMD_marker(data2d)

    Config.loads(data3d)
    st.NMD_marker(data3d)
