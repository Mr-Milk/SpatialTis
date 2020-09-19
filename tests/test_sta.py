import os

import pytest
from anndata import read_h5ad

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.WORKING_ENV = None


def test_read_data(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    pytest.data = data


def test_cell_components():
    data = pytest.data
    pytest.data = data
    st.cell_components(data)
    st.cell_components(data, export_key=CONFIG.cell_components_key, return_df=True)
    CONFIG.WORKING_ENV = "jupyter_notebook"
    st.cell_components(data)
    CONFIG.WORKING_ENV = None


def test_cell_components_plot():
    data = pytest.data
    sp.cell_components(data,
                       ["Patient", "Part"],
                       selected_types=["9", "11", "7", "5"],
                       sort_type="9",
                       group_order={"Patient": ["HPAP005", "HPAP002"]},
                       palette="Category10",
                       title="sth",
                       xaxis_title="sth",
                       yaxis_title="sth",
                       size=(600, 600),
                       display=False,)
    sp.cell_components(data,
                       ["Patient"],
                       selected_types=["9", "11", "7", "5"],
                       group_order={"Patient": ["HPAP005", "HPAP002"]},
                       display=False,)
    sp.cell_components(data,
                       ["Patient", "Part"],
                       ascending=False,
                       direction="horizontal",
                       save="test.png",
                       return_plot=True)
    os.remove('test.png')


@pytest.mark.xfail
def test_bar_plot_failed_param_groupby():
    data = pytest.data
    sp.cell_components(data, ["Patient", "T", "P", "G"],)


@pytest.mark.xfail
def test_bar_plot_failed_param_direction():
    data = pytest.data
    sp.cell_components(data, ["Patient"], direction="what?")


def test_cell_density():
    data = pytest.data
    st.cell_density(data, [(100, 100) for _ in range(18)])
    st.cell_density(data, (100, 100), export_key=CONFIG.cell_density_key, return_df=True)


@pytest.mark.xfail
def test_cell_density_failed():
    data = pytest.data
    st.cell_density(data, ("100", "100"))


def test_cell_density_plot():
    data = pytest.data
    sp.cell_density(data,
                    ["Patient"],
                    group_order={"Patient": ["HPAP005", "HPAP002"]},
                    xaxis_title="sth",
                    yaxis_title="sth",
                    size=(600, 600),
                    title="sth")
    sp.cell_density(data,
                    ["Patient"],
                    direction="horizontal",
                    title="cell density",
                    palette="Category10",
                    save="test.png",
                    return_plot=True,)


@pytest.mark.xfail
def test_violin_plot_failed_param_groupby():
    data = pytest.data
    sp.cell_density(data, ["Patient", "T", "P", "G"],)


@pytest.mark.xfail
def test_violin_plot_failed_param_direction():
    data = pytest.data
    sp.cell_density(data, ["Patient"], direction="what?")


def test_cell_co_occurrence():
    data = pytest.data
    st.cell_co_occurrence(data)
    st.cell_co_occurrence(data, export_key=CONFIG.cell_co_occurrence_key, return_df=True)


def test_cell_co_occurrence_plot():
    data = pytest.data
    sp.cell_co_occurrence(data, use="dot", display=False)
    sp.cell_co_occurrence(data, ["Patient", "Part"], use="heatmap", display=False)


def test_cell_morphology():
    data = pytest.data
    st.cell_morphology(data)


def test_cell_morphology_plot():
    data = pytest.data
    sp.cell_morphology(data, display=False)
