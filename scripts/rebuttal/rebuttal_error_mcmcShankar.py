# This script aims at finding the error that starts from the measured value
# In previous scripts, the mean of the error samples are different from the 
# measured value, bad-looking.
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from skimage import io,morphology,measure,filters
import h5py

organelles = [
    "peroxisome",
    "vacuole",
    "ER",
    "golgi",
    "mitochondria",
    "LD"

] 


# %%
# both functions do not have return values, i.e.:
# their executions are not expected to be assigned to some variables
# instead, the first input will be altered by the functions.
def random_sample(args):
    # return np.random.random(args)
    return 0.1 + 0.8*np.random.random(args)

def upsample(mask,prob):
    dilated = morphology.binary_dilation(mask)
    edge = np.logical_xor(mask,dilated)
    to_compare = prob[edge]
    randoms  = random_sample(to_compare.shape)
    compared = (to_compare > randoms)
    mask[edge] = compared
    return None

def downsample(mask,prob):
    eroded = morphology.binary_erosion(mask)
    edge = np.logical_xor(mask,eroded)
    to_compare = prob[edge] # not (1 - prob[edge]), because last line means a flip
    randoms  = random_sample(to_compare.shape)
    compared = (to_compare > randoms)
    mask[edge] = compared
    return None

# %%
N_sample = 20000
stem = "EYrainbow_glu-100_field-3"

path_cell = f"images/cell/EYrainbow_glucose_largerBF/binCell_{stem}.tif"
img_cell = io.imread(str(path_cell))

dfs = []
for organelle in organelles:
    print(f"Start processing: {organelle}")
    path_organelle = f"images/preprocessed/EYrainbow_glucose_largerBF/probability_{organelle}_{stem}.h5"
    with h5py.File(str(path_organelle)) as h5:
        img_orga = h5["exported_data"][1]

    index = []
    trues = []
    means = []
    stdvs = []
    for cell in measure.regionprops(img_cell):
        min_row, min_col, max_row, max_col = cell.bbox
        img_cell_crop = cell.image
        
        img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col]
        for z in range(img_orga_crop.shape[0]):
            img_orga_crop[z] = img_orga_crop[z] * img_cell_crop

        mask_selected = (img_orga_crop > 0.5)
        mask_dynamic  = np.copy(mask_selected)
        sizes = np.empty(N_sample)
        sizes[0] = np.count_nonzero(mask_selected)
        for i in range(N_sample-1):
            seed = np.random.random()
            if seed < 0.5:
                downsample(mask_dynamic,img_orga_crop)
            else:
                upsample(  mask_dynamic,img_orga_crop)
            sizes[i+1] = np.count_nonzero(mask_dynamic)
    
        index.append(cell.label)
        trues.append(sizes[0])
        means.append(sizes.mean())
        stdvs.append(sizes.std())

        print(f"... simulated cell #{cell.label}")
    dfs.append(pd.DataFrame({
        "organelle": organelle,
        "index"    : index,
        "segmented": trues,
        "average"  : means,
        "standard_deviation": stdvs
    }))
    print(f"Finished: {organelle}")

    
df = pd.concat(dfs,ignore_index=True)
df.to_csv("plots/rebuttal_error/mcmcShankar_10-90.csv",index=False)

# %%
df = pd.read_csv("plots/rebuttal_error/mcmcShankar_25-75.csv")
# %%
df["std/mean"]  = df["standard_deviation"]/df["average"]
df["diff"]      = df["segmented"] - df["average"]
df["diff/mean"] = df["diff"]/df["average"]
df = df[df['diff/mean'].lt(200)]
df.dropna(inplace=True)

# %%
by_organlle = df[["organelle","std/mean","diff/mean"]].groupby("organelle").mean()
plt.figure()
plt.errorbar(
    np.arange(6),np.zeros(6),fmt='None',
    yerr=by_organlle.loc[organelles,'std/mean'],
    capsize=5,ecolor='k'
)
plt.scatter(
    np.arange(6),by_organlle.loc[organelles,"diff/mean"],
    c='k'
)
plt.xticks(
    ticks=np.arange(6),
    labels=organelles
)
plt.xlabel("Organelle")
plt.ylabel("Segmentation Error")
# plt.show()
plt.savefig("plots/rebuttal_error/mcmcShankar_25-75_summary.png",dpi=300)
# %%
for organelle in organelles:
    df_organelle = df[df["organelle"].eq(organelle)]
    plt.figure()
    plt.errorbar(
        np.arange(len(df_organelle)),
        np.zeros(len(df_organelle)),
        yerr=df_organelle["standard_deviation"]/df_organelle["average"],
        capsize=5
    )
    plt.scatter(
        np.arange(len(df_organelle)),
        df_organelle["diff"]/df_organelle["average"]
    )
    plt.title(organelle)
    plt.ylim(-0.5,0.5)
    plt.show()

    print(
        organelle,
        len(df_organelle),
        np.count_nonzero(
            np.absolute(df_organelle["diff"]) > df_organelle["standard_deviation"]
        )
    )


# %% [markdown]
# THE FOLLOWING CELLS CREATES DEMO OF THE ABOVE SIMULATIONS

# %%
img_prob = np.zeros((60,60),dtype=bool)
img_prob[ 7:22, 7:22] = morphology.disk(radius=7, dtype=bool)
img_prob[30:51,30:51] = morphology.disk(radius=10,dtype=bool)
img_prob = filters.gaussian(img_prob,sigma=3)
img_prob = (img_prob - img_prob.min())/(img_prob.max() - img_prob.min())
# %%
img_mask = (img_prob>0.5)

img_external = morphology.binary_dilation(img_mask)
img_external = np.logical_xor(img_mask,img_external)

img_internal = morphology.binary_erosion(img_mask)
img_internal = np.logical_xor(img_mask,img_internal)

# %%
new_down = np.copy(img_mask)
downsample(new_down,img_prob)

new_up = np.copy(img_mask)
upsample(new_up,img_prob)

# %%
