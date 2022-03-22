import os
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData


class FakeSpatialData:
    cell_id = None
    cell_type = None
    gene_name = None
    exp = None
    x = None
    y = None

    cell_meta = None
    cell_loc = None
    exp_matrix = None

    p = None

    def __init__(self, cell_num, gene_num, type_num, save_dir):
        if isinstance(save_dir, str):
            p = Path(save_dir)
        else:
            p = save_dir
        if not p.exists():
            os.mkdir(p)

        self.p = p

        self.cell_id = [f"cell{i}" for i in range(0, cell_num)]

        type_pool = [f'type{i}' for i in range(type_num)]
        self.cell_type = np.random.choice(type_pool, cell_num)

        self.gene_name = [f'gene{i}' for i in range(gene_num)]

        self.x = np.random.randint(0, 10000, cell_num) / 10
        self.y = np.random.randint(0, 10000, cell_num) / 10

        self.cell_loc = pd.DataFrame({'cell_id': self.cell_id, 'x': self.x, 'y': self.y})

        self.cell_meta = pd.DataFrame({'cell_id': self.cell_id,
                                       'cell_type': self.cell_type,
                                       'roi': ['roi' for i in range(cell_num)],
                                       'centroid': [str((x_, y_)) for x_, y_ in zip(self.x, self.y)]})

        self.exp = [np.random.randn(cell_num) for _ in range(gene_num)]
        self.exp_matrix = pd.DataFrame(data=self.exp, columns=self.cell_id)
        self.exp_matrix.insert(0, 'gene', self.gene_name)

    def to_file(self):

        self.cell_meta.to_csv(self.p / "meta.txt", sep="\t", index=False)
        self.cell_loc.to_csv(self.p / "loc.txt", sep="\t", index=False)
        self.exp_matrix.to_csv(self.p / "exp.txt", sep="\t", index=False)

    def to_anndata(self):
        data = AnnData(X=np.asarray(self.exp).T, obs=self.cell_meta, var=pd.DataFrame({"gene": self.gene_name}))
        data.write(self.p / "data.h5ad")


cell_num = [1000]
gene_num = 100
type_num = 30

for c in cell_num:
    try:
        os.mkdir("fake_data")
        os.mkdir("result")
    except Exception:
        pass
    data = FakeSpatialData(c, 100, 20, f'fake_data/data_{c}')
    data.to_file()
    data.to_anndata()
