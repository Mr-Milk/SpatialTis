from anndata import read_h5ad

import spatialtis.plotting as sp
import spatialtis.sta as st
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_sta(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    st.cell_components(data, selected_types=["4", "5", "6", "7"])
    sp.cell_components(data, ["Patient", "Part"])

    st.cell_density(data, (1000, 1000), selected_types=["4", "5", "6", "7"])
    sp.cell_density(data, ["Patient"], direction="horizontal", title="cell density")

    st.cell_co_occurrence(data)
    sp.cell_co_occurrence(data, ["Patient", "Part"])

    st.cell_morphology(data)
    sp.cell_morphology(data, "Patient", "Part", "4")
