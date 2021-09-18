# Since we run in docker, there is not way to detect the actual cpus used by docker
# Had to manually set it
import sys
from time import time

import anndata as ad
import pandas as pd
from memory_profiler import memory_usage

import spatialtis as st

CPUs = int(sys.argv[1])

st.Config.cell_type_key = 'cell_type'
st.Config.centroid_key = 'centroid'
st.Config.exp_obs = ['roi']
st.Config.marker_key = 'gene'
st.Config.verbose = False

cell_num = [1000, 2000, 5000, 10000, 50000]

nns_time = dict(software=[],
                cell_num=[],
                exec_time=[],
                exec_mem=[],
                cpu_count=[],
                )

cci_time = dict(software=[],
                cell_num=[],
                exec_time=[],
                exec_mem=[],
                cpu_count=[],
                )


def time_neighbors(data):
    s1 = time()
    st.find_neighbors(data, expand=3)
    s2 = time()
    return s2 - s1


def time_na(data):
    s1 = time()
    st.cell_interaction(data, order=False, resample=500)
    s2 = time()
    return s2 - s1


heat = False
mem_prof = dict(max_usage=True, multiprocess=True, retval=True)

for c in cell_num:
    data = ad.read_h5ad(f"fake_data/data_{c}/data.h5ad")

    if not heat:
        time_neighbors(data)
        time_na(data)
        heat = True

    for _ in range(3):
        nns_exec_mem, nns_exec_time = memory_usage((time_neighbors, (data,)), **mem_prof)
        cci_exec_mem, cci_exec_time = memory_usage((time_na, (data,)), **mem_prof)

        nns_time['software'].append('spatialtis')
        nns_time['cell_num'].append(c)
        nns_time['exec_time'].append(nns_exec_time)
        nns_time['exec_mem'].append(nns_exec_mem)
        nns_time['cpu_count'].append(CPUs)

        cci_time['software'].append('spatialtis')
        cci_time['cell_num'].append(c)
        cci_time['exec_time'].append(cci_exec_time)
        cci_time['exec_mem'].append(cci_exec_mem)
        cci_time['cpu_count'].append(CPUs)

pd.DataFrame(nns_time).to_csv(f"result/spatialtis_nns_time_{CPUs}core.csv", index=False)
pd.DataFrame(cci_time).to_csv(f"result/spatialtis_cci_time_{CPUs}core.csv", index=False)
