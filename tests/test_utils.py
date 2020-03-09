from anndata import read_h5ad

from spatialtis.utils import prepare_svca

from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"

data = read_h5ad("tests/data/small.h5ad")

def test_svca():
    prepare_svca(data, "tmp/")
