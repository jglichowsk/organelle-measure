# %%
import pandas as pd
import os
from pathlib import Path
from skimage import io,measure
from organelle_measure.tools import batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

def parse_meta_cell(name):
    """name is the stem of the ORGANELLE label image file."""
    #assign various labels according to file naming convention.
    labels=name.split('_')
    experiment=labels[1]
    condition=labels[2]
    field=labels[3]  
    return {
        "experiment": experiment,
        "condition":  condition,
        "hour":       3,
        "field":      field,
    }

def measure1cell(path_in,path_out):
    img_cell = io.imread(str(path_in))
    name = Path(path_in).stem
    meta = parse_meta_cell(name)
    measured = measure.regionprops_table(
                    img_cell,
                    properties=('label','area','centroid','bbox','eccentricity') ### image_intensity returns ndarray of pixels inside bounding box
               )                                                                 ### similarly intensity_mean and intensity_std except those return floats.
    result = meta | measured
    df = pd.DataFrame(result)
    df.rename(columns={'label':'idx-cell'},inplace=True)
    df.to_csv(str(path_out),index=False)
    return None

# %%
imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 

list_in = []
list_out = []

for folder in folders:
    if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/cell_measure')):
        print('Creating folder.')
        os.makedirs(newpath)
    else:
        print(str(folder+'/cell_measure'),'already there.')

    for path_in in (Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment'))).glob("binCell*.tif"):
        path_out = Path(str(imgs+'/'+exp+'/'+folder+'/cell_measure'))/f"{path_in.stem.partition('-')[2]}.csv"
        list_in.append(path_in)
        list_out.append(path_out)


args = pd.DataFrame({
    "path_in":   list_in,
    "path_out":  list_out
})

batch_apply(measure1cell,args)

# %%

# %%
# import numpy as np
# import matplotlib.pyplot as plt

# px_x,px_y,px_z = 0.41,0.41,0.20
# csv["effective-volume"] = (px_x*px_y)*np.sqrt(px_x*px_y)*(2.)*csv.loc[:,"area"]*np.sqrt(csv.loc[:,"area"])/np.sqrt(np.pi)

# plt.figure()
# plt.hist(csv["area"],bins=20)
# plt.xlabel("size/pixels")
# plt.savefig("binCell-after_EYrainbow_glu-0-5_field-0_histogram-pixel.png")

# plt.figure()
# plt.hist(csv["effective-volume"],bins=20)
# plt.xlabel("size/microns")
# plt.savefig("binCell-after_EYrainbow_glu-0-5_field-0_histogram-microns.png")
