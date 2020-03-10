from anndata import read_h5ad
from pathlib import Path
from spatialtis import CONFIG
from spatialtis.utils import prepare_svca
DATA_DIR = str(Path("tests/data/small.h5ad").absolute())

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"

data = read_h5ad(DATA_DIR)


def test_svca():
    prepare_svca(data, "tmp/")
