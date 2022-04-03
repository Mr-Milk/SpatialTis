import pytest

import spatialtis as st
from spatialtis import Config


# Spatial analysis
def test_cell_dispersion(data_2d):
    st.cell_dispersion(data_2d, r=50, method="id")
    st.cell_dispersion(data_2d, quad=(10, 10), method="morisita")
    st.cell_dispersion(data_2d, method="clark_evans")


def test_cell_dispersion_3d(data_3d):
    st.cell_dispersion(data_3d, r=50, method="id")


def test_spatial_hetero(data_2d):
    st.spatial_heterogeneity(data_2d, method="shannon")
    st.spatial_heterogeneity(data_2d, method="leibovici", d=10)
    st.spatial_heterogeneity(data_2d, method="altieri", cut=3)


def test_spatial_hetero_3d(data_3d):
    st.spatial_heterogeneity(data_3d, method="shannon")
    st.spatial_heterogeneity(data_3d, method="leibovici", d=10)
    st.spatial_heterogeneity(data_3d, method="altieri", cut=3)


def test_hotspot(data_2d):
    st.hotspot(data_2d, rect_side=(10, 10))


# find neighbors
def test_find_neighbors(data_2d):
    st.find_neighbors(data_2d, method="delaunay")
    st.find_neighbors(data_2d, r=30, k=5, method="kdtree")
    # st.find_neighbors(data_2d, scale=1.2, method="rtree")


def test_find_neighbors_3d(data_3d):
    st.find_neighbors(data_3d, r=30, k=5, method="kdtree")


def test_spatial_autocorr(data_2d):
    st.spatial_autocorr(data_2d, method="moran_i")
    st.spatial_autocorr(data_2d, method="geary_c")


def test_spatial_autocorr_3d(data_3d):
    st.spatial_autocorr(data_3d, method="moran_i")
    st.spatial_autocorr(data_3d, method="geary_c")





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
