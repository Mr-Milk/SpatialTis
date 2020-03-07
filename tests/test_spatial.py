from anndata import read_h5ad

import spatialtis.plotting as sp
import spatialtis.spatial as ss
from spatialtis import CONFIG, Neighbors

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"

data = read_h5ad("../tmp/small.h5ad")

ss.spatial_distribution(data, r=50)
ss.spatial_distribution(data, quad=(10, 10), method="MID", overwrite=True)
ss.spatial_distribution(data, method="CE", overwrite=True)

ss.spatial_heterogeneity(data, compare=0)
ss.hotspot(data, grid_size=10, selected_types=["4", "6"])

n = Neighbors(data)
n.find_neighbors(expand=5)
n.neighbors_count()
n.export_neighbors()
n.read_neighbors()

ss.communities(n)
ss.neighborhood_analysis(n, resample=50)
ss.spatial_enrichment_analysis(n, resample=50)

sp.spatial_distribution(data, ["Patient", "Part"])
sp.spatial_heterogeneity(data, ["Patient"])
