import os
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
import anndata as ad
import spatialtis as st
from spatialtis import Config

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


@pytest.fixture(scope="session")
def ome_images():
    images_dir = TEST_DIR / "ome_images"
    images = [
        images_dir / "slide0_ROI7.ome.tiff",
        images_dir / "slide0_ROI8.ome.tiff"
    ]
    masks = [
        images_dir / "slide0_ROI7_mask.ome.tif",
        images_dir / "slide0_ROI8_mask.ome.tif"
    ]
    return images, masks


@pytest.fixture(scope="session")
def visium_folder():
    return [
        TEST_DIR / 'visium_data' / 'GSM4644079_v9_filtered_gene_bc_matrix_spatial',
        TEST_DIR / 'visium_data' / 'GSM4644080_v10_filtered_gene_bc_matrix_spatial',
    ]


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
    c = Config
    c.reset()
    c.roi_key = "ROI"
    c.cell_type_key = "cell_type"
    if (loc == 'obs') & (dtype == 'wkt'):
        data.obs[['x', 'y']] = points
        st.wkt_points(data, ['x', 'y'], write_config=False)

        c.centroid_key = "centroid"
        c.dumps(data)
    elif loc == 'obsm':
        data.obsm['spatial'] = points
        c.dumps(data)
    else:
        if is3D:
            data.obs[['x', 'y', 'z']] = points
            c.centroid_key = ['x', 'y', 'z']
            c.dumps(data)
        else:
            data.obs[['x', 'y']] = points
            c.centroid_key = ['x', 'y']
            c.dumps(data)
    return data


@pytest.fixture(scope="session")
def data_shape():
    data = ad.read_h5ad(TEST_DIR / "tmp" / "small.h5ad")
    c = Config
    c.reset()
    c.roi_key = 'ROI_ID'
    c.shape_key = 'cell_shape'
    c.cell_type_key = 'leiden'
    c.marker_key = 'Markers'
    Config.dumps(data)
    return data


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

