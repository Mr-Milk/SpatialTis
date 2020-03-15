import pandas as pd
import numpy as np

from anndata import read_h5ad

import spatialtis.plotting as sp
from spatialtis import CONFIG

CONFIG.EXP_OBS = ["Patients", "Sample", "ROI"]
CONFIG.CELL_TYPE_COL = "leiden"

# synthetic data
mindex = pd.MultiIndex.from_arrays([
    np.hstack((['Ted' for t in range(0, 20)], ['Gia' for t in range(0, 20)])),
    np.hstack((['2015' for t in range(0, 10)], ['2016' for t in range(0, 10)], ['2015' for t in range(0, 10)],
               ['2016' for t in range(0, 10)])),
    ['Store' if t % 2 == 0 else 'Online' for t in range(0, 40)],
    [i for i in range(0, 40)]],
    names=('Manager', 'Year', 'Sales', 'id'))

col = ('apple', 'banana', 'cherry', 'blueberry', 'pineapple', 'strawberry')

df = pd.DataFrame(columns=col, index=mindex)
for label, _ in df.items():
    df[label] = 40 * np.random.randn(40) + 200
df.columns.names = ['fruits']


def test_cell_map(shared_datadir):
    data = read_h5ad(shared_datadir / 'small.h5ad')
    ROI = data.obs[
        (data.obs["Patient"] == "HPAP005")
        & (data.obs["Part"] == "Tail")
        & (data.obs["ROI"] == "ROI1")
        ]
    sp.cell_map(ROI)


def test_violin():
    sp.violin_plot(df, ['Manager', 'Year', 'Sales'], 'apple', split='Sales')


def test_bar():
    sp.stacked_bar(df, ['Manager', 'Year'])
