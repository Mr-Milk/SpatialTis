from anndata import read_h5ad
from spatialtis import CONFIG
from spatialtis.utils import prepare_svca

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_svca(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    prepare_svca(data, "tmp/")
