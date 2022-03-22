# Since we run in docker, there is not way to detect the actual cpus used by docker
# Had to manually set it
import sys
from time import time

import anndata as ad
import numpy as np
import pandas as pd
from memory_profiler import memory_usage

import squidpy as sq

CPUs = int(sys.argv[1])
heat = False
K = 10
TIMES = 1000

dataset = ('Stereo-seq-MouseEmbryo', 'IMC-BreastCancer')
mature = ad.read_h5ad("data/E16-5.h5ad")
basel = ad.read_h5ad("data/IMC-scp-basel.h5ad")

nns_tb = dict(software=[],
              dataset=[],
              cpu_count=[],
              exec_time=[],
              exec_mem=[],
              )

cci_tb = dict(software=[],
              dataset=[],
              cpu_count=[],
              exec_time=[],
              exec_mem=[],
              )


def time_neighbors(data):
    s1 = time()
    sq.gr.spatial_neighbors(data, n_neighs=K)
    s2 = time()
    print(f"Find neighbors used {s2-s1:.2}s")
    return s2 - s1


def time_na(data):
    s1 = time()
    sq.gr.nhood_enrichment(data, "annotation", n_perms=TIMES)
    s2 = time()
    print(f"Cell interaction used {s2 - s1:.2}s")
    return s2 - s1


def time_autocorr(data):
    s1 = time()
    sq.gr.spatial_autocorr(data, mode="moran", two_tailed=True)
    s2 = time()
    return s2 - s1


def time_all(data):
    s1 = time()
    sq.gr.spatial_neighbors(data, n_neighs=K)
    sq.gr.nhood_enrichment(data, "cell_type", n_perms=TIMES, show_progress_bar=False)
    s2 = time()
    return s2 - s1


mem_prof = dict(max_usage=True, multiprocess=True, retval=True)

# Test Stereo-seq data
nns_mem, nns_time = memory_usage((time_neighbors, (mature,)), **mem_prof)
print("Finished finding neighbor")
cci_mem, cci_time = memory_usage((time_na, (mature,)), **mem_prof)
print("Finished cell interaction")
at_mem, at_time = memory_usage((time_autocorr, (mature,)), **mem_prof)
print(f"Finished spatial autocorr, {at_mem} {at_time}")

nns_tb['software'].append('Squidpy')
nns_tb['dataset'].append(dataset[0])
nns_tb['exec_time'].append(nns_time)
nns_tb['exec_mem'].append(nns_mem)
nns_tb['cpu_count'].append(CPUs)

cci_tb['software'].append('Squidpy')
cci_tb['dataset'].append(dataset[0])
cci_tb['exec_time'].append(cci_time)
cci_tb['exec_mem'].append(cci_mem)
cci_tb['cpu_count'].append(CPUs)

pd.DataFrame(nns_tb).to_csv(f"result/squidpy_embryo_nns_{CPUs}core.csv", index=False)
pd.DataFrame(cci_tb).to_csv(f"result/squidpy_embryo_cci_{CPUs}core.csv", index=False)

# Test IMC dataset, Only run when CPUS >= 8
# if CPUs >= 8:
#     multi_mem = []
#     multi_time = []
#     failed = 0
#     for c in basel.obs['core'].unique():
#         roi = basel[basel.obs['core'] == c]
#         try:
#             roi_mem, roi_time = memory_usage((time_all, (roi,)), **mem_prof)
#             multi_mem.append(roi_mem)
#             multi_time.append(roi_time)
#         except:
#             failed += 1
#     print(f"Processing Failed: {failed}")
#     pd.DataFrame({
#         "software": ["Squidpy"],
#         "dataset": [dataset[1]],
#         "cpu_count": [CPUs],
#         "exec_time": [np.sum(multi_time)],
#         "exec_mem": [np.amax(multi_mem)]
#     }).to_csv(f"result/squidpy_multi_all_{CPUs}cores.csv", index=False)
