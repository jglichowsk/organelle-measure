# This script aims at finding the error that starts from the measured value
# In previous scripts, the mean of the error samples are different from the 
# measured value, bad-looking.
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from skimage import io,morphology,measure,filters,segmentation
import sys
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
# from pathing_variables import expmt_path
from tools import batch_apply

organelles = [
    # "px",
    # "vo",
    # "er",
    # "gl",
    # "mt",
    "ld"
] 

org_path = r'_______' #path to folder containing all organelle images to be analyzed
def parse_meta_organelle(name: str):
    """
    Args: Name is the stem of the ORGANELLE label image file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention.
    deconv, stk, time, organelle, date, strain, condition, field=name.split('_')
    # field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
        "time":       time[-1]
    }
# %%
# both functions do not have return values, i.e.:
# their executions are not expected to be assigned to some variables
# instead, the first input will be altered by the functions.
def random_sample(args: tuple):
    """
    Args: Tuple representing shape of current fov/crop
    Outputs: Array of same shape with values between 0.25 and 0.75
    """
    # return np.random.random(args)
    return 0.25 + 0.5*np.random.random(args)

def upsample(mask: np.ndarray, prob: np.ndarray):
    """
    Args: Organelle mask and probability arrays
    Outputs: None, but expands/erodes org mask
    """
    dilated = morphology.binary_dilation(mask)
    edge = np.logical_xor(mask,dilated)
    to_compare = prob[edge]
    randoms  = random_sample(to_compare.shape)
    compared = (to_compare > randoms)
    mask[edge] = compared
    return None

def downsample(mask: np.ndarray, prob: np.ndarray):
    """
    Args: Organelle mask and probability arrays
    Outputs: None, but expands/erodes org mask
    """
    eroded = morphology.binary_erosion(mask)
    edge = np.logical_xor(mask,eroded)
    to_compare = prob[edge] # not (1 - prob[edge]), because last line means a flip
    randoms  = random_sample(to_compare.shape)
    compared = (to_compare > randoms)
    mask[edge] = compared
    return None

# %%
def simulate1fov(folder_in: str, folder_out: str, stem: str, N_sample: int=1):
    """
    Args: folder_in = folder containing all expmt subfolders; folder_out = where to save output csv; 
            stem = identifier to distinguish different runs?; N_sample = # sim runs
    Outputs: None, but saves csv containing mcmc analysis metrics
    """
    # N_sample = 20000 #Number of simulation runs

    # finished = []
    # finished = [p.stem.partition("25-75_")[2] for p in Path("plots/rebuttal_error/optimization").glob("*.csv")]
    # if stem in finished:
    #     print(f"Already finished: {stem}")
    #     return
    
    # segm_files = [p for p in Path(folder_in+'/cell_segment').glob("*_segm.tif")]
    # stems = [p.stem[:-9] for p in segm_files]

    ##########FIX THIS#############
    conf_interval=45 #time interval between confocal captures

    meta = parse_meta_organelle(stem)
    path_cell = Path(folder_in+'/cell_segment')/f"BF-timelapse_{stem.split('_', 4)[-1]}_segm-afftransf-expand.tif"
    img_cell = io.imread(path_cell)
    img_cell = img_cell.astype(np.int64)
    # for j in range(len(segm_files)): #for each fov in folder
    #     img_cell = io.imread(segm_files[j])

    dfs = []

    for organelle in organelles: #for each organelle in fov
        print("Start processing: " f"{meta['date']}_frame-{meta['time']}_{organelle}")
        # path_organelle = Path(folder_in)/f"ilastik_prob/{organelle}_{stem}_Probabilities.tiff"
        path_organelle = Path(folder_in)/f"ilastik_prob/{stem}_Probabilities.tiff"
        img_orga = io.imread(path_organelle)

        # path_reference= Path(folder_in)/f"preprocess/{organelle}_{stem}.tif"
        path_reference= Path(folder_in)/f"preprocess/{stem}.tif"
        img_reference  = io.imread(path_reference)

        index = []
        total_trues = []
        total_means = []
        total_stdvs = []

        count_trues = []
        count_means = []
        count_stdvs = []

        mean_trues  = []
        mean_means  = []
        mean_stdvs  = []

        for cell in measure.regionprops(img_cell[(int(meta['time'])-1)*conf_interval, :, :]):
            min_row, min_col, max_row, max_col = cell.bbox
            # print(min_row, min_col, max_row, max_col)
            img_cell_crop = cell.image
            # print(img_cell_crop.shape)
            img_orga_crop = img_orga[:, min_row:max_row, min_col:max_col]
            # print(img_orga.shape)
            # print(img_orga_crop.shape)
            for z in range(img_orga_crop.shape[0]):
                img_orga_crop[z] = img_orga_crop[z] * img_cell_crop #take only org pixels within cell segm mask

            mask_selected = (img_orga_crop > 0.5) #initially, pixels with >0.5 probability counted as true
            mask_dynamic  = np.copy(mask_selected)

            samples_total = np.empty(N_sample)
            samples_total[0] = np.count_nonzero(mask_selected)
            samples_count = np.empty(N_sample)

            if organelle in ["px","gl","ld"]: #for some reason only do this for globular signals
                img_ref_crop = img_reference[:,min_row:max_row,min_col:max_col]
                mask_selected = segmentation.watershed(-img_ref_crop,mask=mask_selected)
            label_organelles = measure.label(mask_selected)

            samples_count[0] = len(measure.regionprops(label_organelles))
            samples_mean = np.empty(N_sample)
            samples_mean[0] = 0 if samples_count[0]==0 else samples_total[0] / samples_count[0]

            for i in range(N_sample-1): #for each mcmc run
                seed = np.random.random()
                if seed < 0.5: #randomly choose to expand or erode org masks
                    downsample(mask_dynamic,img_orga_crop)
                else:
                    upsample(mask_dynamic,img_orga_crop)
                samples_total[i+1] = np.count_nonzero(mask_dynamic)

                # watershed, count, and average
                if organelle in ["px","gl","ld"]: #again w this
                    mask_dynamic = segmentation.watershed(-img_ref_crop,mask=mask_dynamic)
                label_organelles = measure.label(mask_dynamic)
                samples_count[i+1] = len(measure.regionprops(label_organelles))
                samples_mean[i+1] = 0 if samples_count[i+1]==0 else samples_total[i+1] / samples_count[i+1]

            index.append(cell.label)
            total_trues.append(samples_total[0])
            total_means.append(samples_total.mean())
            total_stdvs.append(samples_total.std())

            count_trues.append(samples_count[0])
            count_means.append(samples_count.mean())
            count_stdvs.append(samples_count.std())

            mean_trues.append(samples_mean[0])
            mean_means.append(samples_mean.mean())
            mean_stdvs.append(samples_mean.std())

            print(f"... simulated cell #{cell.label}")
        dfs.append(pd.DataFrame({
            "organelle": organelle,
            "index"    : index,
            "total_true":  total_trues,
            "total_mean":  total_means,
            "total_stdv":  total_stdvs,
            "count_trues": count_trues,
            "count_means": count_means,
            "count_stdvs": count_stdvs,
            "mean_trues":  mean_trues,
            "mean_means":  mean_means,
            "mean_stdvs":  mean_stdvs,
        }))
        print(f"Finished: {organelle}")

    df = pd.concat(dfs,ignore_index=True)
    df.to_csv(f"{folder_out}/mcmc_25-75_{stem}.csv",index=False)
    
    return None

# %%
folder_in=r'C:\Users\jglic\Downloads\07152026 3c-ey2795'
stems = []
for file in Path(folder_in+'/preprocess').glob("*.tif"):
    meta=parse_meta_organelle(file.stem)
    if meta['organelle'] in organelles:
        stems.append(file.stem)
args = pd.DataFrame({
    "folder_in":  [folder_in],
    "folder_out": [folder_in],
    "stem": [stems],
})
# %%
batch_apply(simulate1fov,args)

# %%
df = pd.read_csv("plots/rebuttal_error/mcmcShankar_25-75.csv")
# %%
df["total_std/mean"]  = df["total_stdv"] / df["total_mean"]
df["total_diff"]      = df["total_true"] - df["total_mean"]
df["total_diff/mean"] = df["total_diff"] / df["total_mean"]

df["count_std/mean"]  = df["count_stdvs"] / df["count_means"]
df["count_diff"]      = df["count_trues"] - df["count_means"]
df["count_diff/mean"] = df["count_diff"]  / df["count_means"]

df["mean_std/mean"]  = df["mean_stdvs"] / df["mean_means"]
df["mean_diff"]      = df["mean_trues"] - df["mean_means"]
df["mean_diff/mean"] = df["mean_diff"]  / df["mean_means"]

# df = df[df['diff/mean'].lt(200)]
df.dropna(inplace=True)

# %%
plt.rcParams["figure.autolayout"]=True
plt.rcParams['font.size'] = '14'
by_organlle = df[["organelle","total_std/mean","total_diff/mean"]].groupby("organelle").mean()
plt.figure()
plt.bar(
    np.arange(3), by_organlle.loc[organelles,'total_std/mean'],
    color='1', edgecolor='k', label="upper error"
)
plt.bar(
    np.arange(3), -by_organlle.loc[organelles,'total_std/mean'],
    color='1', edgecolor='k', hatch='/', label="lower error"
)
# plt.scatter(
#     np.arange(3),by_organlle.loc[organelles,"total_diff/mean"],
#     c='k',marker="x", label="measured value"
# )
plt.legend(fontsize=14)
plt.xticks(
    ticks=np.arange(3),
    labels=organelles,
    fontsize=11
)
plt.title("Total Volume")
plt.xlabel("Organelle",fontsize='20')
plt.ylabel("Segmentation Error",fontsize='20')
# plt.show()
plt.savefig("plots/rebuttal_error/3organelles_mcmcMore_25-75_bar_total.png",dpi=600)

# %%
plt.rcParams["figure.autolayout"]=True
plt.rcParams['font.size'] = '14'
by_organlle = df[["organelle","count_std/mean","count_diff/mean"]].groupby("organelle").mean()
plt.figure()
plt.bar(
    np.arange(3), by_organlle.loc[organelles,'count_std/mean'],
    color='1', edgecolor='k', label="upper error"
)
plt.bar(
    np.arange(3), -by_organlle.loc[organelles,'count_std/mean'],
    color='1', edgecolor='k', hatch='/', label="lower error"
)
# plt.scatter(
#     np.arange(3),by_organlle.loc[organelles,"count_diff/mean"],
#     c='k',marker="x", label="measured value"
# )
plt.legend(fontsize=14)
plt.xticks(
    ticks=np.arange(3),
    labels=organelles,
    fontsize=11
)
plt.title("Organelle Population")
plt.xlabel("Organelle",fontsize='20')
plt.ylabel("Segmentation Error",fontsize='20')
# plt.show()
plt.savefig("plots/rebuttal_error/3organelles_mcmcMore_25-75_bar_count.png",dpi=600)

# %%
plt.rcParams["figure.autolayout"]=True
plt.rcParams['font.size'] = '14'
by_organlle = df[["organelle","mean_std/mean","mean_diff/mean"]].groupby("organelle").mean()
plt.figure()
plt.bar(
    np.arange(3), by_organlle.loc[organelles,'mean_std/mean'],
    color='1', edgecolor='k', label="upper error"
)
plt.bar(
    np.arange(3), -by_organlle.loc[organelles,'mean_std/mean'],
    color='1', edgecolor='k', hatch='/', label="lower error"
)
# plt.scatter(
#     np.arange(3),by_organlle.loc[organelles,"mean_diff/mean"],
#     c='k',marker="x", label="measured value"
# )
plt.legend(fontsize=14)
plt.xticks(
    ticks=np.arange(3),
    labels=organelles,
    fontsize=11
)
plt.title("Average Volume")
plt.xlabel("Organelle",fontsize='20')
plt.ylabel("Segmentation Error",fontsize='20')
# plt.show()
plt.savefig("plots/rebuttal_error/3organelles_mcmcMore_25-75_bar_average.png",dpi=600)
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

# %% have a look at the variance of organelle volume fraction across same group of cells
from organelle_measure.data import read_results
px_x,px_y,px_z = 0.41,0.41,0.20

df_bycell = read_results(Path("./data"),["EYrainbow_glucose_largerBF"],(px_x,px_y,px_z))

# %%
df_measure = df_bycell.loc[
                df_bycell["condition"].eq(100)
              & df_bycell["field"].eq(3)
            ]
# %%
total_std = {}
for organelle in organelles:
    total_std[organelle] = df_measure.loc[df_measure["organelle"].eq(organelle),"total"].std()/(px_x * px_y * px_z)


# %%
