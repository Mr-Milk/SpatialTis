import pandas as pd
import pytest

from spatialtis import Config
from spatialtis.utils import col2adata_obs, df2adata_uns
from spatialtis.utils.log import pretty_time


class FakeType:
    a = 0


def test_read_config_property():
    Config
    Config.verbose
    Config.roi_key
    Config.exp_obs
    Config.cell_type_key
    Config.env


def test_read_set_roi_key():
    Config.exp_obs = ["Patient", "Part", "ROI"]
    assert Config.roi_key == "ROI"
    Config.roi_key = "Part"
    assert Config.exp_obs == ["Patient", "ROI", "Part"]


@pytest.mark.xfail
def test_set_roi_key_failed():
    Config.exp_obs = None
    Config.exp_obs = ["Patient", "Part", "ROI"]
    Config.roi_key = "ttt"


def test_set_exp_obs():
    Config.exp_obs = 1
    Config.exp_obs = 10.1
    Config.exp_obs = None
    Config.exp_obs = "ttt"


@pytest.mark.xfail
def test_set_exp_obs_failed():
    ft = FakeType
    Config.exp_obs = ft


def test_set_cell_type_key():
    Config.cell_type_key = 1
    Config.cell_type_key = 10.1
    Config.cell_type_key = None
    Config.cell_type_key = "type"


@pytest.mark.xfail
def test_set_cell_type_key_failed():
    Config.cell_type_key = [123]


def test_set_working_env():
    Config.env = "what"
    Config.env = "jupyter"


def test_set_verbose():
    Config.verbose = True
    Config.verbose = False


@pytest.mark.xfail
def test_set_verbose_failed():
    Config.verbose = 1


def test_pretty_time():
    assert pretty_time(0.01) == "10ms"
    assert pretty_time(1) == "1s0ms"
    assert pretty_time(1.1) == "1s100ms"
    assert pretty_time(61) == "1m1s"
    assert pretty_time(3661) == "1h1m1s"


Config.exp_obs = ["Patient", "Part", "ROI"]
Config.cell_type_key = "leiden"
Config.MARKER_KEY = "Markers"
Config.env = None


def test_data2adata(data):
    df = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 5, 6],
            "b": [4, 5, 6, 7, 8, 9, 2],
        }
    )

    col = [0 for _ in range(len(data.obs))]
    Config.verbose = True
    df2adata_uns(df, data, "test_df")
    col2adata_obs(col, data, "test_col")
    Config.env = "jupyter_notebook"
    df2adata_uns(df, data, "test_df")
    col2adata_obs(col, data, "test_col")
