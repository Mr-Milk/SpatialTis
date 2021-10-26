import pytest

import spatialtis as st
from spatialtis import Config


# Basic analysis
def test_basic(data):
    Config.exp_obs = ["Patient", "Part", "ROI"]
    Config.cell_type_key = "leiden"
    Config.centroid_key = "centroid"
    Config.marker_key = "Markers"
    Config.shape_key = "cell_shape"

    st.transform_points(data, 'centroid')
    st.transform_shapes(data, 'cell_shape')

    st.cell_components(data)
    st.cell_density(data)
    st.cell_co_occurrence(data)
    st.cell_morphology(data)


# Spatial analysis
def test_cell_dispersion(data):
    st.cell_dispersion(data, r=50, method="id")
    st.cell_dispersion(data, quad=(10, 10), method="morisita")
    st.cell_dispersion(data, method="clark_evans")


def test_spatial_hetero(data):
    st.spatial_heterogeneity(data, method="shannon")
    st.spatial_heterogeneity(data, method="leibovici")
    st.spatial_heterogeneity(data, method="altieri")


def test_hotspot(data):
    st.hotspot(data, rect_side=(10, 10))


# find neighbors
def test_find_neighbors(data):
    st.find_neighbors(data, method="delaunay")
    st.find_neighbors(data, r=30, k=5, method="kdtree")
    st.find_neighbors(data, scale=1.2, method="rtree")


def test_spatial_autocorr(data):
    st.spatial_autocorr(data, method="moran_i")
    st.spatial_autocorr(data, method="geary_c")


def test_somde(data):
    st.somde(data, epoch=10)


def test_neighborhood(data):
    st.cell_interaction(data)
    st.cell_interaction(data, method="zscore")
    st.cell_interaction(data, order=True, export_key="neighborhood_order")


def test_spatial_enrichment(data):
    st.spatial_enrichment(data)
    st.spatial_enrichment(data, threshold=0.1)
    st.spatial_enrichment(data, order=True, export_key="enrichment_order")


def test_spatial_coexp(data):
    st.spatial_coexp(data, method="pearson")
    st.spatial_coexp(data, method="spearman")


def test_ncd_markers(data):
    st.NCD_marker(data)


def test_nmd_markers(data):
    st.NMD_marker(data)


def test_svca(data, tmpdir):
    st.prepare_svca(data, tmpdir)


@pytest.mark.xfail
def test_svca_failed(data, tmpdir):
    st.prepare_svca(data, tmpdir, marker_key="wrong_key")
