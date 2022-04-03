import pandas as pd
import pytest

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config


ROI = "1"
TYPES = ["12", "14", "8"]
MARKERS = ['CD20', 'CD57', 'C-peptide', 'Ghrelin']
GORDER = {"Patient": ["HPAP005", "HPAP002"]}


def test_set_config(tmpdir):
    Config.exp_obs = ["Patient", "Part", "ROI_ID"]
    Config.cell_type_key = "leiden"
    Config.marker_key = "Markers"
    Config.shape_key = "cell_shape"
    Config.centroid_key = "centroid"
    Config.env = None
    Config.auto_save = True
    Config.save_path = tmpdir


def test_cell_map(data):
    sp.cell_map(data, ROI, use_shape=True)
    sp.cell_map(data, ROI, selected_types=TYPES)


def test_expression_map(data):
    sp.expression_map(data, ROI, "CD20")


def test_neighbors_map(data):
    sp.neighbors_map(data, ROI)


def test_cell_components(data):
    sp.cell_components(data, "Patient")


def test_cell_density(data):
    sp.cell_density(data, "Patient")


def test_cell_morphology(data):
    sp.cell_morphology(data, 'Patient')


def test_cell_co_occurrence(data):
    sp.cell_co_occurrence(data)


def test_cell_dispersion(data):
    sp.cell_dispersion(data)


def test_cell_interaction(data):
    sp.cell_interaction(data, order=True)
    sp.cell_interaction(data)


def test_spatial_enrichment(data):
    sp.spatial_enrichment(data, order=True)
    sp.spatial_enrichment(data)

