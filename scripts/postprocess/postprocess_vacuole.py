# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,morphology,measure
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import skeletonize_zbyz,watershed_zbyz,find_complete_rings,better_vacuole_img,batch_apply


# https://stackoverflow.com/questions/28281742/fitting-a-circle-to-a-binary-image


def postproc_vacuole(path_in: str,path_cell: str,path_out: str, threshold=0.5):
    bkgdprob=io.imread(path_in)
    orgprob=1-bkgdprob
    img_org = (orgprob>threshold)
    # img_org = np.argmax(orgprob,axis=0)
    # img_org=(img_org>0)
    # img_maxslice=np.argmax(org_prob,axis=0)

    img_skeleton = skeletonize_zbyz(img_org)
    img_core = find_complete_rings(img_skeleton)
    img_cell = io.imread(str(path_cell))[0,:,:]
    if (img_core.shape[1]+img_core.shape[2]) != (img_cell.shape[0]+img_cell.shape[1]): #fix disparity between camera and confocal detector img sizing
        img_mask = np.zeros((img_core.shape[1],img_core.shape[2]),dtype=int) 
        shape0,shape1 = img_cell.shape
        img_mask[:shape0,:shape1] = img_cell
        img_cell=img_mask
    
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
list_in=[]
list_cell=[] 
list_out=[]


if not os.path.exists(newpath:=Path(expmt_path+'/postprocess')):
    print('Creating folder ',str(expmt_path+'/postprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/ilastik_prob').glob("vo*.tiff"): 
    path_parts=path_in.stem.split("_")
    path_end="_".join(path_parts[:-1])
    cell_parts=path_in.stem.split('-')
    cell_end="-".join(cell_parts[:3])[3:]

    path_out=Path(newpath)/f'{path_end}_label.tif'
    path_ref=Path(expmt_path+'/preprocess')/f"{path_end}.tif"
    path_cell=Path(expmt_path+'/cell_segment')/f"BF-timelapse_{cell_end}_segm.tif"

    list_in.append(path_in)
    list_cell.append(path_cell)
    list_out.append(path_out)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})
# %%
batch_apply(postproc_vacuole,args)

# %%
