 from anndata import read_h5ad

import spatialtis.plotting as sp
import spatialtis.spatial as ss
from spatialtis import CONFIG, Neighbors

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_spatial_dist(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ss.spatial_distribution(data, r=50)
    ss.spatial_distribution(data, quad=(10, 10), method="quad")
    ss.spatial_distribution(data, method="nns")

    sp.spatial_distribution(data, ["Patient", "Part"])


def test_spatial_hetero(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ss.spatial_heterogeneity(data, compare=0)

    sp.spatial_heterogeneity(data, ["Patient"])


def test_hotspot(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ss.hotspot(data, grid_size=10)


def test_neighbors(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    n = Neighbors(data)
    n.find_neighbors(expand=5)
    n.neighbors_count()
    n.export_neighbors()
    n.read_neighbors()


def test_neighbor_analysis(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    n = Neighbors(data)
    n.find_neighbors(expand=5)
    ss.communities(n)
    ss.neighborhood_analysis(n, resample=50)
    ss.spatial_enrichment_analysis(n, resample=50)

    sp.neighborhood_analysis(data, ["Patient", "Part"])
    sp.spatial_enrichment_analysis(data, ["Patient", "Part"])


def test_spatial_mp(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    n = Neighbors(data)
    n.find_neighbors(scale=1.1, mp=True)
    ss.neighborhood_analysis(n, resample=50, mp=True)
    ss.spatial_enrichment_analysis(n, resample=50, mp=True)
    ss.hotspot(data, grid_size=10)
