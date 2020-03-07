from anndata import read_h5ad

from spatialtis.utils import prepare_svca

data = read_h5ad("../tmp/small.h5ad")

prepare_svca(data, "../tmp/")
