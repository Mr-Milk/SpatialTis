import os
import shutil
from pathlib import Path

import matplotlib
import pytest
from anndata import read_h5ad

matplotlib.use('agg')

TEST_DIR = Path(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session", autouse=True)
def create_data():
    test_data_dir = TEST_DIR / "data"
    tmp_dir = TEST_DIR / "tmp"
    shutil.copytree(test_data_dir, tmp_dir)
    yield
    shutil.rmtree(tmp_dir)


@pytest.fixture(scope="module")
def tmpdir():
    return TEST_DIR / "tmp"


@pytest.fixture(scope="session")
def data():
    return read_h5ad(TEST_DIR / "tmp" / "small.h5ad")
