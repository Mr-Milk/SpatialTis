import pytest
from anndata import read_h5ad

from spatialtis import CONFIG
from spatialtis.utils import prepare_svca, pretty_time, timer


@timer(prefix="say sth before", suffix="say sth after")
def fake_func(groupby=None):
    pass


class FakeType:
    a = 0


def test_read_config_property():
    CONFIG
    CONFIG.VERBOSE
    print(CONFIG.VERBOSE)
    CONFIG.ROI_KEY
    CONFIG.EXP_OBS
    CONFIG.CELL_TYPE_KEY
    CONFIG.WORKING_ENV
    CONFIG.CPU_ALLOC


def test_read_set_roi_key():
    CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
    CONFIG.ROI_KEY = "ROI"
    CONFIG.ROI_KEY = "Part"


@pytest.mark.xfail
def test_set_roi_key_failed():
    CONFIG.EXP_OBS = None
    CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
    CONFIG.ROI_KEY = "ttt"


def test_set_exp_obs():
    CONFIG.EXP_OBS = 1
    CONFIG.EXP_OBS = 10.1
    CONFIG.EXP_OBS = None
    CONFIG.EXP_OBS = "ttt"


@pytest.mark.xfail
def test_set_exp_obs_failed():
    ft = FakeType
    CONFIG.EXP_OBS = ft


def test_set_cell_type_key():
    CONFIG.CELL_TYPE_KEY = 1
    CONFIG.CELL_TYPE_KEY = 10.1
    CONFIG.CELL_TYPE_KEY = None
    CONFIG.CELL_TYPE_KEY = "type"


@pytest.mark.xfail
def test_set_cell_type_key_failed():
    CONFIG.CELL_TYPE_KEY = [123]


def test_set_working_env():
    CONFIG.WORKING_ENV = "what"
    CONFIG.WORKING_ENV = "jupyter_notebook"
    CONFIG.WORKING_ENV = "jupyter_lab"
    CONFIG.WORKING_ENV = "nteract"
    CONFIG.WORKING_ENV = "zeppelin"


def test_set_cpu_alloc():
    CONFIG.CPU_ALLOC = 1


@pytest.mark.xfail
def test_set_cpu_alloc_failed():
    CONFIG.CPU_ALLOC = 1.2


def test_set_verbose():
    CONFIG.VERBOSE = True
    CONFIG.VERBOSE = False


@pytest.mark.xfail
def test_set_verbose_failed():
    CONFIG.VERBOSE = 1


def test_pretty_time():
    assert pretty_time(0.01) == "10ms"
    assert pretty_time(1) == "1s0ms"
    assert pretty_time(1.1) == "1s100ms"
    assert pretty_time(61) == "1m1s"
    assert pretty_time(3661) == "1h1m1s"


@pytest.mark.xfail
def test_timer_failed():
    fake_func()


CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.MARKER_KEY = "Markers"
CONFIG.WORKING_ENV = None


def test_timer():
    fake_func(groupby="what")


def test_svca(shared_datadir, tmpdir):
    CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
    CONFIG.CELL_TYPE_KEY = "leiden"
    data = read_h5ad(shared_datadir / 'small.h5ad')
    prepare_svca(data, tmpdir)
