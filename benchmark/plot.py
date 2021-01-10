import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")

st_nns_1 = pd.read_csv("spatialtis_nns_time_1core.csv")
go_nns_1 = pd.read_csv("giotto_nns_time_1core.csv")
st_nns_4 = pd.read_csv("spatialtis_nns_time_4core.csv")
go_nns_4 = pd.read_csv("giotto_nns_time_4core.csv")
st_nns_1['core'] = "1 core"
go_nns_1['core'] = "1 core"
st_nns_4['core'] = "4 cores"
go_nns_4['core'] = "4 cores"
nns_time = pd.concat([st_nns_1, go_nns_1, st_nns_4, go_nns_4])
nns_time['software'] = nns_time['software'].str.cat(nns_time['core'], sep="-")
nns_time['cell_num'] = nns_time['cell_num'].astype('str')

nns = sns.barplot(x="cell_num", y="exec_time", hue="software", data=nns_time, ci=95)
nns.set(title="Neighbor search execute time", ylabel="Time (s)", xlabel="Number of cell")
nns.get_figure().savefig("compare_neighbor_search.png", dpi=300, bbox_inches="tight")

st_cci_1 = pd.read_csv("spatialtis_cci_time_1core.csv")
go_cci_1 = pd.read_csv("giotto_cci_time_1core.csv")
st_cci_4 = pd.read_csv("spatialtis_cci_time_4core.csv")
go_cci_4 = pd.read_csv("giotto_cci_time_4core.csv")
st_cci_1['core'] = "1 core"
go_cci_1['core'] = "1 core"
st_cci_4['core'] = "4 cores"
go_cci_4['core'] = "4 cores"
cci_time = pd.concat([st_cci_1, go_cci_1, st_cci_4, go_cci_4])
cci_time['software'] = cci_time['software'].str.cat(cci_time['core'], sep="-")
cci_time['cell_num'] = cci_time['cell_num'].astype('str')

cci = sns.barplot(x="cell_num", y="exec_time", hue="software", data=cci_time, ci=95)
cci.set(title="Cell-cell interaction execute time", ylabel="Time (s)", xlabel="Number of cell")
cci.get_figure().savefig("compare_cell_cell_interaction.png", dpi=300, bbox_inches="tight")
