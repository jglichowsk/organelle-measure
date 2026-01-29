import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from organelle_measure.data import read_results

px_x,px_y,px_z = 0.41,0.41,0.20
organelles = [
    "peroxisome",
    "vacuole",
    "ER",
    "golgi",
    "mitochondria",
    "LD"
]


subfolders = [
    "EYrainbow_glucose_largerBF",
]
bycell_6color = read_results(Path("./data/results"),subfolders,(px_x,px_y,px_z))
bycell_6color = bycell_6color.loc[
									bycell_6color["condition"].eq(100)
								 ]
bycell_6color.rename(
	columns={
		"idx_cell"      : "cell-idx",
		"cell-volume"   : "cell-volume-um3",
		"mean"          : "mean-um3",
		"total"         : "total-um3",
		"total-fraction": "volume-fraction"
	},
	inplace=True
)


# 1-color cells
dfs_read = []
for path in Path("data/results/2024-02-16_rebuttal1color").glob("*.csv"):
	dfs_read.append(pd.read_csv(str(path)))
df_1color = pd.concat(dfs_read,ignore_index=True)
df_1color["volume-um3"] = (px_x*px_y*px_z)*df_1color["area"]
df_1color.loc[df_1color["organelle"].eq("vacuole"),"volume-um3"] = (
	(px_x * px_y * df_1color.loc[df_1color["organelle"].eq("vacuole"),"area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1color.loc[df_1color["organelle"].eq("vacuole"),"area"]
		/np.pi
	)
)
df_1color["cell-volume-um3"] = (
	(px_x * px_y * df_1color["cell-area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1color["cell-area"]
		/np.pi
	)
)

bycell_1color = df_1color[["organelle","cell-idx","cell-volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1color["mean-um3"]  = df_1color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1color["total-um3"] = df_1color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).sum()
bycell_1color["count"]     = df_1color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).count()
bycell_1color["volume-fraction"] = bycell_1color["total-um3"] / bycell_1color["cell-volume-um3"]
bycell_1color.reset_index(inplace=True)


# 3-color cells
dfs_read = []
for path in Path("data/results/2024-02-06_paperRebuttal3Colors").glob("*.csv"):
	dfs_read.append(pd.read_csv(str(path)))
df_3color = pd.concat(dfs_read,ignore_index=True)
df_3color["volume-um3"] = (px_x*px_y*px_z)*df_3color["area"]
df_3color.loc[df_3color["organelle"].eq("vacuole"),"volume-um3"] = (
	(px_x * px_y * df_3color.loc[df_3color["organelle"].eq("vacuole"),"area"])
	*2
	*np.sqrt(
		px_x * px_y * df_3color.loc[df_3color["organelle"].eq("vacuole"),"area"]
		/np.pi
	)
)
df_3color["cell-volume-um3"] = (
	(px_x * px_y * df_3color["cell-area"])
	*2
	*np.sqrt(
		px_x * px_y * df_3color["cell-area"]
		/np.pi
	)
)

bycell_3color = df_3color[["organelle","cell-idx","cell-volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_3color["mean-um3"]  = df_3color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_3color["total-um3"] = df_3color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).sum()
bycell_3color["count"]     = df_3color[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).count()
bycell_3color["volume-fraction"] = bycell_3color["total-um3"] / bycell_3color["cell-volume-um3"]
bycell_3color.reset_index(inplace=True)



# Plot
for organelle in organelles:
	for property in ["mean-um3","total-um3","volume-fraction","count"]:
		values = []
		errors = []
		for dataset in [bycell_6color,bycell_3color,bycell_1color]:
			values.append(dataset.loc[dataset["organelle"].eq(organelle),property].mean())
			errors.append(dataset.loc[dataset["organelle"].eq(organelle),property].std()/np.sqrt(len(dataset.loc[dataset["organelle"].eq(organelle),property])))
		plt.figure()
		plt.bar(
			np.arange(3), height=values, yerr=errors,

		)
		plt.xticks(
			ticks=np.arange(3),
			labels=["All 6-label diploids","3-label diploids","1-color diploids"]
		)
		name = property.replace('um3','volume').replace("count","number").replace('-',' ').title()
		y_label = r"Volume ($\mu$m$^3$)" if "um" in property else name
		plt.ylabel(y_label)
		plt.title(f"{organelle}\n{name}")
		plt.savefig(
			f"data/compare_diploids/{organelle}_{name}.png",
			dpi=600
		)
		

