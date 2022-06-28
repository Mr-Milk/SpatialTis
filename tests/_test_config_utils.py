import pandas as pd
import pytest

from spatialtis import Config
from spatialtis.utils import col2adata, df2adata_uns
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
    Config.reset()
    Config.exp_obs = ["Patient", "Part", "ROI_ID"]
    Config.roi_key = "Part"
    assert Config.exp_obs == ["Patient", "ROI_ID", "Part"]
    Config.exp_obs = ["Patient", "Part", "ROI_ID"]
    assert Config.roi_key == "ROI_ID"


def test_config_dumps_loads(data_shape):
    Config.dumps(data_shape)
    Config.loads(data_shape)
    Config.reset()


@pytest.mark.xfail
def test_set_roi_key_failed():
    Config.reset()
    Config.exp_obs = ["Patient", "Part", "ROI"]
    Config.roi_key = "ttt"


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


Config.exp_obs = ["Patient", "Part", "ROI_ID"]
Config.cell_type_key = "leiden"
Config.MARKER_KEY = "Markers"
Config.env = None

