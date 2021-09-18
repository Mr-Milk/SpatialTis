import pytest

import spatialtis as st
from spatialtis import Config


# Basic analysis
def test_basic(data):
    Config.exp_obs = ["Patient", "Part", "ROI"]
    Config.cell_type_key = "leiden"
    Config.MARKER_KEY = "Markers"
    Config.MP = False

    st.cell_components(data)
    st.cell_density(data)
    st.cell_co_occurrence(data)
    st.cell_morphology(data)


# Spatial analysis
def test_spatial_dist(data):
    st.spatial_distribution(
        data,
        r=50,
        method="vmr",
    )
    st.spatial_distribution(data, quad=(10, 10), method="quad")
    st.spatial_distribution(data, method="nns")


def test_spatial_dist_mp(data):
    st.spatial_distribution(data, r=50, method="vmr", mp=True)
    st.spatial_distribution(data, quad=(10, 10), method="quad", mp=True)
    st.spatial_distribution(data, method="nns", mp=True)


def test_spatial_hetero(data):
    st.spatial_heterogeneity(data, method="shannon", compare="Patient")
    st.spatial_heterogeneity(data, method="leibovici")
    st.spatial_heterogeneity(data, method="altieri")


def test_spatial_hetero_mp(data):
    st.spatial_heterogeneity(data, method="shannon", compare="Patient", mp=True)
    st.spatial_heterogeneity(data, method="leibovici", mp=True)
    st.spatial_heterogeneity(data, method="altieri", mp=True)


def test_hotspot(data):
    st.hotspot(data, grid_size=10)


def test_hotspot_mp(data):
    st.hotspot(data, grid_size=10, mp=True)


# find neighbors
def test_find_neighbors(data):
    st.find_neighbors(data, expand=10)
    st.find_neighbors(data, expand=10, use_shape=True)


@pytest.mark.xfail
def test_find_neighbors_failed(data):
    st.find_neighbors(data, scale=2)


def test_neighborhood(data):
    st.cell_interaction(data)
    st.cell_interaction(data, method="zscore")
    st.cell_interaction(data, order=True, export_key="neighborhood_order")


def test_spatial_enrichment(data):
    st.spatial_enrichment(data)
    st.spatial_enrichment(data, threshold=0.1)
    st.spatial_enrichment(data, order=True, export_key="enrichment_order")


def test_spatial_coexpression(data):
    st.spatial_coexp(data, method="pearson")
    st.spatial_coexp(
        data, selected_markers=data.var[Config.MARKER_KEY].tolist()
    )


def test_community(data):
    st.cell_community(data)


def test_ncd_markers(data):
    st.NCDMarkers(data)


def test_nmd_markers(data):
    st.NMDMarkers(data)


def test_svca(data, tmpdir):
    st.prepare_svca(data, tmpdir)


@pytest.mark.xfail
def test_svca_failed(data, tmpdir):
    st.prepare_svca(data, tmpdir, marker_key="wrong_key")
