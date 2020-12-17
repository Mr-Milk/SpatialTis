import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")

st_nns = pd.read_csv("spatialtis_nns_time.csv")
go_nns = pd.read_csv("giotto_nns_time.csv")
nns_time = pd.concat([st_nns, go_nns])

nns = sns.barplot(x="cell_num", y="exec_time", hue="software", data=nns_time, ci=95).set_title("Neighbor search execute time")
nns.get_figure().savefig("nns.png")

st_cci = pd.read_csv("spatialtis_cci_time.csv")
go_cci = pd.read_csv("giotto_cci_time.csv")
cci_time = pd.concat([st_cci, go_cci])

cci = sns.barplot(x="cell_num", y="exec_time", hue="software", data=cci_time, ci=95).set_title("Cell-cell interaction execute time")
cci.get_figure().savefig("cci.png")
