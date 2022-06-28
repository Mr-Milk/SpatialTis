import pytest

import spatialtis as st
from spatialtis import Config


# Spatial analysis
@pytest.mark.parametrize("method,kw", [
    ("id", dict(r=50)),
    ("morisita", dict(quad=(10, 10))),
    ("clark_evans", dict())
])
def test_cell_dispersion(data_2d, method, kw):
    st.cell_dispersion(data_2d, method=method, **kw)


def test_cell_dispersion_3d(data_3d):
    st.cell_dispersion(data_3d, r=50, method="id")


@pytest.mark.parametrize("method,kw", [
    ("leibovici", dict(d=10)),
    ("altieri", dict(cut=3)),
    ("shannon", dict())
])
def test_spatial_hetero(data_2d, data_3d, method, kw):
    st.spatial_heterogeneity(data_2d, method=method, **kw)
    st.spatial_heterogeneity(data_3d, method=method, **kw)


def test_hotspot(data_2d):
    st.hotspot(data_2d, rect_side=(10, 10))


# find neighbors
@pytest.mark.parametrize("method,kw", [
    ("kdtree", dict(r=30)),
    ("kdtree", dict(k=5)),
    ("kdtree", dict(k=5, r=30)),
    ("delaunay", dict())
])
def test_find_neighbors(data_2d, method, kw):
    st.find_neighbors(data_2d, method=method, **kw)
    # st.find_neighbors(data_2d, scale=1.2, method="rtree")


@pytest.mark.parametrize("method,kw", [
    ("kdtree", dict(r=30)),
    ("kdtree", dict(k=5)),
    ("kdtree", dict(k=5, r=30)),
])
def test_find_neighbors_3d(data_3d, method, kw):
    st.find_neighbors(data_3d, method=method, **kw)


@pytest.mark.parametrize("method", ["moran_i", "geary_c"])
def test_spatial_autocorr(data_2d, data_3d, method):
    st.spatial_autocorr(data_2d, method=method)
    st.spatial_autocorr(data_3d, method=method)


def test_neighborhood(data_2d):
    st.cell_interaction(data_2d)
    st.cell_interaction(data_2d, method="zscore")
    st.cell_interaction(data_2d, export_key="neighborhood_order")


def test_spatial_enrichment(data_2d):
    st.spatial_enrichment(data_2d)
    st.spatial_enrichment(data_2d, threshold=0.1)
    st.spatial_enrichment(data_2d, export_key="enrichment_order")


def test_spatial_coexp(data_2d):
    st.spatial_coexp(data_2d)
    st.spatial_coexp(data_2d, use_cell_type=True)


def test_ncd_markers(data_2d):
    st.NCD_marker(data_2d)


def test_nmd_markers(data_2d):
    st.NMD_marker(data_2d)
