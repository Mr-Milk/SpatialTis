import os
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import anndata as ad
import spatialtis as st

matplotlib.use('agg')

TEST_DIR = Path(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session", autouse=True)
def create_data():
    test_data_dir = TEST_DIR / "data"
    tmp_dir = TEST_DIR / "tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    shutil.copytree(test_data_dir, tmp_dir)
    yield
    shutil.rmtree(tmp_dir)


@pytest.fixture(scope="session")
def tmpdir():
    return TEST_DIR / "tmp"


def fake(roi, N, markers, is3D=False, loc='obs', dtype="wkt"):

    dims = 3 if is3D else 2
    points = np.random.randint(1, 100, (N * roi, dims))
    roi_names = []
    for i, _ in enumerate(range(roi)):
        roi_names += [f"ROI_{i + 1}" for _ in range(N)]

    markers_name = ["".join(np.random.choice(list("abcdefghijklmn"), 3)) for _ in range(markers)]

    obs = pd.DataFrame({"ROI": roi_names, "cell_type": np.random.choice(list("abcdefg"), N * roi)})
    var = pd.DataFrame(index=markers_name)
    X = np.random.randn(N * roi, markers)
    data = ad.AnnData(obs=obs, var=var, X=X)
    if (loc == 'obs') & (dtype == 'wkt'):
        data.obs[['x', 'y']] = points
        st.transform_points(data, ['x', 'y'], write_config=False)
    elif loc == 'obsm':
        data.obsm['spatial'] = points
    else:
        if is3D:
            data.obs[['x', 'y', 'z']] = points
        else:
            data.obs[['x', 'y']] = points
    return data


@pytest.fixture(scope="session")
def data_shape():
    return ad.read_h5ad(TEST_DIR / "tmp" / "small.h5ad")


@pytest.fixture(scope="session")
def data_wkt():
    return fake(roi=2, N=100, markers=10, is3D=False)


@pytest.fixture(scope="session")
def data2d_keys():
    return fake(roi=2, N=100, markers=10, is3D=False, loc="obs", dtype='int')


@pytest.fixture(scope="session")
def data3d_keys():
    return fake(roi=2, N=100, markers=10, is3D=True, loc="obs", dtype='int')


@pytest.fixture(scope="session")
def data2d():
    return fake(roi=2, N=100, markers=10, is3D=False, loc="obsm")


@pytest.fixture(scope="session")
def data3d():
    return fake(roi=2, N=100, markers=10, is3D=True, loc="obsm")

