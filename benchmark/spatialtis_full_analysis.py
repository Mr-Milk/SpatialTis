from time import time

import anndata as ad

import spatialtis as st
from spatialtis import CONFIG

data = ad.read_h5ad("../data/imc_data.h5ad")
CONFIG.EXP_OBS = ["stage", "case", "part", "image"]
CONFIG.CELL_TYPE_KEY = "cell_type"
CONFIG.MARKER_KEY = "markers"
CONFIG.CENTROID_KEY = "centroid"

t1 = time()
# Basic Analysis
st.cell_components(data)
st.cell_density(data)
st.cell_co_occurrence(data)
st.cell_morphology(data, metric_key="eccentricity")

# Spatial Analysis
islets_cells = ['gamma', 'delta', 'alpha', 'beta']
st.spatial_distribution(data)
st.hotspot(data)
st.spatial_heterogeneity(data)

st.find_neighbors(data, expand=8)
st.cell_community(data)
st.neighborhood_analysis(data)

# exclude DNA markers
selected_markers = ["INS", "CD38", "CD44", "PCSK2", "CD99", "CD68", "MPO", "SLC2A1",
                    "CD20", "AMY2A", "CD3e", "PPY", "PIN", "GCG", "PDX1", "SST", "SYP", "KRT19",
                    "CD45", "FOXP3", "CD45RA", "CD8a", "IAPP", "NKX6-1", "CD4", "PTPRN", "cCASP3"]

st.spatial_co_expression(data, selected_markers=selected_markers)
st.NCDMarkers(data, use_cell_type=True, selected_markers=selected_markers)

t2 = time()
print(f"Running full analysis of SpatialTis on IMC dataset with {data.shape} using {round(t2 - t1, 3)}s")
