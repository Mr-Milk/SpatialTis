from time import time

import anndata as ad
import pandas as pd

import spatialtis as st

st.CONFIG.CELL_TYPE_KEY = 'cell_type'
st.CONFIG.CENTROID_KEY = 'centroid'
st.CONFIG.EXP_OBS = ['roi']
st.CONFIG.MARKER_KEY = 'gene'
st.CONFIG.VERBOSE = False

cell_num = [1000, 5000, 10000]#, 50000, 100000]

nns_time = dict(software=[],
                cell_num=[],
                exec_time=[])

cci_time = dict(software=[],
                cell_num=[],
                exec_time=[])


def time_neighbors(n):
    s1 = time()
    n.find_neighbors(expand=3)
    s2 = time()
    return s2 - s1


def time_na(n):
    s1 = time()
    st.neighborhood_analysis(n, order=False)
    s2 = time()
    return s2 - s1


heat = False

for c in cell_num:
    data = ad.read_h5ad(f"data_{c}/data.h5ad")
    n = st.Neighbors(data, 'point')

    if not heat:
        time_neighbors(n)
        time_na(n)
        heat = True

    for _ in range(3):
        nns_exec_time = time_neighbors(n)
        cci_exec_time = time_na(n)

        nns_time['software'].append('spatialtis')
        nns_time['cell_num'].append(c)
        nns_time['exec_time'].append(nns_exec_time)

        cci_time['software'].append('spatialtis')
        cci_time['cell_num'].append(c)
        cci_time['exec_time'].append(cci_exec_time)

pd.DataFrame(nns_time).to_csv("spatialtis_nns_time.csv", index=False)
pd.DataFrame(cci_time).to_csv("spatialtis_cci_time.csv", index=False)
