# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,morphology,measure
from organelle_measure.tools import skeletonize_zbyz,watershed_zbyz,find_complete_rings,better_vacuole_img,batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

def postproc_vacuole(path_in,path_cell,path_out):
    bkgdprob=io.imread(str(path_in))
    orgprob=1-bkgdprob
    # img_orga = np.argmax(org_prob,axis=0)
    # img_maxslice=np.argmax(org_prob,axis=0)
    img_org = (orgprob>bkgdprob)

    img_cell = io.imread(str(path_cell))
    img_skeleton  = skeletonize_zbyz(img_org)

    img_core      = find_complete_rings(img_skeleton)
    
    # img_vacuole   = better_vacuole_img(img_core,img_watershed)
    img_vacuole = np.zeros_like(img_core,dtype=int)
    for z in range(img_vacuole.shape[0]):
        sample = img_core[z]
        candidates = np.unique(sample[img_cell>0])
        for color in candidates:
            if len(np.unique(img_cell[sample==color]))==1:
                img_vacuole[z,sample==color] = color

    io.imsave(
        str(path_out),
        util.img_as_uint(img_vacuole) 
    )
    return None

# %%
list_in=[]; list_cell=[]; list_out=[]; #path lists for function batch process
imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on.

for folder in folders:
    if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/postprocess')):
        print('Creating folder.')
        os.makedirs(newpath)
    else:
        print(str(folder+'/postprocess'),'already there.')

    for path_cell in (Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment'))).glob(f"*.tif"):
        path_binary = (Path(str(imgs+'/'+exp+'/'+folder+'/ilastik prob')))/f"vacuole_{path_cell.stem.partition('-')[2]}_Probabilities.tiff"
        path_output = (Path(str(newpath)))/f"label-vacuole_{path_cell.stem.partition('-')[2]}.tiff"
        list_in.append(path_binary)
        list_cell.append(path_cell)
        list_out.append(path_output)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})
# args.to_csv("./vauocle.csv",index=False)

batch_apply(postproc_vacuole,args)

# %%
