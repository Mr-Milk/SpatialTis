from anndata import read_h5ad

from spatialtis import CONFIG
from spatialtis.utils import prepare_svca

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_KEY = "leiden"
CONFIG.WORKING_ENV = None


def test_svca(shared_datadir, tmpdir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    prepare_svca(data, tmpdir)
