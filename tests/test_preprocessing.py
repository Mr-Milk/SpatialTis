from spatialtis.preprocessing import read_ROIs
import pytest

DATA_ENTRY = "data/shrink_data/Patients"
METADATA = "data/shrink_data/metadata.csv"

conditions = ['Patients', 'Sample', 'ROI']

data = read_ROIs(DATA_ENTRY, conditions)
data.config_file(METADATA, channel_col='channels', marker_col='markers')

adata = data.to_anndata(mp=True)
