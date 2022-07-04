from __future__ import annotations

import h5py
import pandas as pd
import warnings
from anndata import AnnData
from pathlib import Path
from scipy.io import mmread
from scipy.sparse import csr_matrix
from typing import List, Optional

from spatialtis.config import _Config


def read_10x_h5(h5file):
    with h5py.File(h5file) as h:
        # csr_matrix((data, indices, indptr), [shape=(M, N)])
        mtx = h.get('matrix')
        shape = mtx.get('shape')[:]
        exp = csr_matrix((mtx.get('data'), mtx.get('indices'), mtx.get('indptr')), shape=(shape[1], shape[0]))

        features = mtx.get('features')
        gene_names = features['name'][:].astype(str)
        gene_ids = features['id'][:].astype(str)
        barcodes = mtx.get('barcodes')[:].astype(str)

        return barcodes, gene_names, gene_ids, exp


def read_10x_folder(path, stfile):
    coord = read_10x_spatial(stfile)
    path = Path(path)
    bc = [i for i in path.glob("*barcodes*")][0]
    features = [i for i in path.glob("*features*")][0]
    mtx = [i for i in path.glob("*matrix.mtx*")][0]

    try:
        bc_data = pd.read_csv(bc, sep="\t", header=None)
    except OSError:
        bc_data = pd.read_csv(bc, sep="\t", header=None, compression=None)

    try:
        f_data = pd.read_csv(features, sep="\t", header=None)
    except OSError:
        f_data = pd.read_csv(features, sep="\t", header=None, compression=None)

    try:
        exp = mmread(mtx)
    except OSError:
        exp = mmread(open(mtx))
    barcodes = bc_data.to_numpy().flatten()
    gene_names = f_data.iloc[:, 1].to_numpy().flatten()
    gene_ids = f_data.iloc[:, 0].to_numpy().flatten()
    # sort the coord based on barcode
    sort_bc = pd.concat([coord, bc_data.set_index(0)], join="inner", axis=1).index
    coord = coord.loc[sort_bc, :]
    if len(sort_bc) != len(barcodes):
        warnings.warn(f"{len(barcodes) - len(sort_bc)} barcodes are missing")
        exp = pd.DataFrame.sparse.from_spmatrix(exp, columns=barcodes).T.loc[sort_bc, :].sparse.to_coo()
    else:
        exp = exp.T
    return sort_bc, gene_names, gene_ids, exp, coord


def read_10x_spatial(stlist):
    coord = pd.read_csv(stlist, header=None, index_col=0)
    coord = coord[coord.iloc[:, 0] == 1]
    return coord.iloc[:, [1, 2]].copy()


def read_visium(
        paths: List[str | Path],
        read_filtered: bool = True,  # "filtered" or "raw"
        annotations: Optional[pd.DataFrame] = None,
) -> AnnData:
    """Read visium data from visium result folders.

    Parameters
    ----------
    paths : list of str or path
            Visium result folders.
    read_filtered : bool, default: True
        To use filtered matrix or raw matrix.
    annotations : pd.DataFrame, default: None
        Your annotations on the ROI, for example:
        if you have two ROI, you can annotate it with `pd.DataFrame({'ROI': ['ROI1', 'ROI2'])`.

    """
    if isinstance(paths, (Path, str)):
        paths = [paths]

    if annotations is None:
        annotations = pd.DataFrame({"ROI": [f"ROI{i + 1}" for i in range(len(paths))]})
    roi_header = annotations.columns

    barcodes = []
    coords = []
    exp_dfs = []
    anno_collect = []
    path_records = []

    for p, (_, names) in zip(paths, annotations.iterrows()):
        p = Path(p)
        h5file = p / "filtered_feature_bc_matrix.h5" if read_filtered else p / "raw_feature_bc_matrix.h5"
        spatial_file = p / 'spatial' / 'tissue_positions_list.csv'
        if not h5file.exists():
            folder = p / "filtered_feature_bc_matrix" if read_filtered else p / "raw_feature_bc_matrix"
            bc, gene_names, gene_ids, exp, coord = read_10x_folder(folder, spatial_file)
        else:
            bc, gene_names, gene_ids, exp = read_10x_h5(h5file)
            coord = read_10x_spatial(spatial_file)

        exp_index = pd.DataFrame({"name": gene_names, "id": gene_ids})
        sparse_exp = pd.DataFrame.sparse.from_spmatrix(data=exp, columns=pd.MultiIndex.from_frame(exp_index))

        coords.append(coord)
        exp_dfs.append(sparse_exp)
        barcodes += bc.tolist()
        anno_collect += [names.values for _ in range(len(coord))]
        path_records += [str(p) for _ in range(len(coord))]

    obs = pd.concat(coords)
    obs.columns = ['array_x', 'array_y']
    obs['barcodes'] = barcodes
    obs['path'] = path_records
    obs = obs.reset_index(drop=True)
    obs.index = obs.index.astype(str)
    obs[roi_header] = anno_collect

    exp_all = pd.concat(exp_dfs).fillna(0)
    var = exp_all.columns.to_frame()
    var.set_index("id")

    data = AnnData(obs=obs,
                   var=var,
                   X=exp_all.sparse.to_coo().tocsr()
                   )

    config = _Config()
    config.exp_obs = annotations.columns.tolist()
    config.centroid_key = 'spatial'
    config.marker_key = 'name'
    config.dumps(data)

    return data
