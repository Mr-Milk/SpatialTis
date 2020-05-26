from anndata import read_h5ad

import spatialtis.plotting as sp
import spatialtis.sta as st
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patient", "Part", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"


def test_sta(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    st.cell_components(data)
    sp.cell_components(data, ["Patient", "Part"])
    sp.cell_components(data, ["Patient", "Part"], direction="horizontal",)

    st.cell_density(data, (100, 100))
    sp.cell_density(data, ["Patient"], title="cell density")
    sp.cell_density(data, ["Patient"], direction="horizontal", title="cell density")

    st.cell_co_occurrence(data)
    sp.cell_co_occurrence(data, ["Patient", "Part"])

    st.cell_morphology(data)
    sp.cell_morphology(data, "Patient", "Part", "4")
