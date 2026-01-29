import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from organelle_measure.data import read_results

px_x,px_y,px_z = 0.41,0.41,0.20
organelles = [
    "peroxisome",
    # "vacuole",
    "ER",
    "golgi",
    "mitochondria",
    "LD"
]

# 6-color cells segmented by 6-color ilastik
subfolders = [
    "EYrainbow_glucose_largerBF"
]
bycell_6segmented_by6 = read_results(Path("./data/results"),subfolders,(px_x,px_y,px_z))
bycell_6segmented_by6 = bycell_6segmented_by6.loc[bycell_6segmented_by6["condition"].eq(100)]
bycell_6segmented_by6.rename(
	columns={
		"idx_cell"      : "cell-idx",
		"cell-volume"   : "cell-volume-um3",
		"mean"          : "mean-um3",
		"total"         : "total-um3",
		"total-fraction": "volume-fraction"
	},
	inplace=True
)


# 1-color cells segmented by 1-color ilastik
dfs_read = []
for path in Path("data/results/2024-02-16_rebuttal1color").glob("*.csv"):
	dfs_read.append(pd.read_csv(str(path)))
df_1segmented_by1 = pd.concat(dfs_read,ignore_index=True)
df_1segmented_by1["volume-um3"] = (px_x*px_y*px_z)*df_1segmented_by1["area"]
df_1segmented_by1.loc[df_1segmented_by1["organelle"].eq("vacuole"),"volume-um3"] = (
	(px_x * px_y * df_1segmented_by1.loc[df_1segmented_by1["organelle"].eq("vacuole"),"area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1segmented_by1.loc[df_1segmented_by1["organelle"].eq("vacuole"),"area"]
		/np.pi
	)
)
df_1segmented_by1["cell-volume-um3"] = (
	(px_x * px_y * df_1segmented_by1["cell-area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1segmented_by1["cell-area"]
		/np.pi
	)
)

bycell_1segmented_by1 = df_1segmented_by1[["organelle","cell-idx","cell-volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1segmented_by1["mean-um3"]  = df_1segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1segmented_by1["total-um3"] = df_1segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).sum()
bycell_1segmented_by1["count"]     = df_1segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).count()
bycell_1segmented_by1["volume-fraction"] = bycell_1segmented_by1["total-um3"] / bycell_1segmented_by1["cell-volume-um3"]
bycell_1segmented_by1.reset_index(inplace=True)


# 6-color cells segmented by 1-color ilastik
dfs_read = []
for path in Path("data/results/2024-02-21_rebuttal1n6color/ilastik-1_color-6").glob("*.csv"):
	df = pd.read_csv(str(path))
	if len(df) > 0:
		dfs_read.append(df)
df_6segmented_by1 = pd.concat(dfs_read,ignore_index=True)
df_6segmented_by1["cell-idx"] = df_6segmented_by1["field"].astype(str) + "-" + df_6segmented_by1["cell-idx"].astype(str)
df_6segmented_by1["volume-um3"] = (px_x*px_y*px_z)*df_6segmented_by1["area"]
df_6segmented_by1.loc[df_6segmented_by1["organelle"].eq("vacuole"),"volume-um3"] = (
	(px_x * px_y * df_6segmented_by1.loc[df_6segmented_by1["organelle"].eq("vacuole"),"area"])
	*2
	*np.sqrt(
		px_x *  px_y * df_6segmented_by1.loc[df_6segmented_by1["organelle"].eq("vacuole"),"area"]
		/np.pi
	)
)
df_6segmented_by1["cell-volume-um3"] = (
	(px_x * px_y * df_6segmented_by1["cell-area"])
	*2
	*np.sqrt(
		px_x * px_y * df_6segmented_by1["cell-area"]
	    /np.pi
	)
)

bycell_6segmented_by1 = df_6segmented_by1[["organelle","cell-idx","cell-volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_6segmented_by1["mean-um3"]  = df_6segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_6segmented_by1["total-um3"] = df_6segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).sum()
bycell_6segmented_by1["count"]     = df_6segmented_by1[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).count()
bycell_6segmented_by1["volume-fraction"] = bycell_6segmented_by1["total-um3"] / bycell_6segmented_by1["cell-volume-um3"]
bycell_6segmented_by1.reset_index(inplace=True)


# 1-color cells segmented by 6-color ilastik
dfs_read = []
for path in Path("data/results/2024-02-21_rebuttal1n6color/ilastik-6_color-1").glob("*.csv"):
	df = pd.read_csv(str(path))
	if len(df) > 0:
		dfs_read.append(df)
df_1segmented_by6 = pd.concat(dfs_read,ignore_index=True)
df_1segmented_by6["volume-um3"] = (px_x*px_y*px_z)*df_1segmented_by6["area"]
df_1segmented_by6.loc[df_1segmented_by6["organelle"].eq("vacuole"),"volume-um3"] = (
	(px_x * px_y * df_1segmented_by6.loc[df_1segmented_by6["organelle"].eq("vacuole"),"area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1segmented_by6.loc[df_1segmented_by6["organelle"].eq("vacuole"),"area"]
		/np.pi
	)
)
df_1segmented_by6["cell-volume-um3"] = (
	(px_x * px_y * df_1segmented_by6["cell-area"])
	*2
	*np.sqrt(
		px_x * px_y * df_1segmented_by6["cell-area"]
	    /np.pi
	)
)

bycell_1segmented_by6 = df_1segmented_by6[["organelle","cell-idx","cell-volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1segmented_by6["mean-um3"]  = df_1segmented_by6[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).mean()
bycell_1segmented_by6["total-um3"] = df_1segmented_by6[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).sum()
bycell_1segmented_by6["count"]     = df_1segmented_by6[["organelle","cell-idx","volume-um3"]].groupby(["organelle","cell-idx"]).count()
bycell_1segmented_by6["volume-fraction"] = bycell_1segmented_by6["total-um3"] / bycell_1segmented_by6["cell-volume-um3"]
bycell_1segmented_by6.reset_index(inplace=True)


# Plot
for organelle in organelles:
	for property in ["mean-um3","total-um3","volume-fraction","count"]:
		plt.figure()
		for d,dataset in enumerate([bycell_6segmented_by6,bycell_6segmented_by1,bycell_1segmented_by1,bycell_1segmented_by6]):
			plt.bar(
			    [d], height=[dataset.loc[dataset["organelle"].eq(organelle),property].mean()], 
			           yerr=[dataset.loc[dataset["organelle"].eq(organelle),property].std()],
			)
		plt.xticks(
			ticks=np.arange(4),
			labels=[
				"6-color cells\n6-color ilastik",
				"6-color cells\n1-color ilastik",
				"1-color cells\n1-color ilastik",
				"1-color cells\n6-color ilastik",
			]
		)
		name = property.replace('um3','volume').replace("count","number").replace('-',' ').title()
		y_label = r"Volume ($\mu$m$^3$)" if "-um" in property else name
		plt.ylabel(y_label)
		plt.title(f"{organelle}\n{name}")
		plt.savefig(
			f"data/compare_diploids/crossvalidate1n6_{organelle}_{name}.png",
			dpi=600
		)
		

