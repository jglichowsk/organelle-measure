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

def postproc_vacuole(path_in: str,path_cell: str,path_out: str, threshold=0.5):
    bkgdprob=io.imread(path_in)
    orgprob=1-bkgdprob
    # img_orga = np.argmax(org_prob,axis=0)
    # img_maxslice=np.argmax(org_prob,axis=0)
    img_org = (orgprob>threshold)

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
list_in=[]
list_cell=[] 
list_out=[]


if not os.path.exists(newpath:=Path(expmt_path+'/postprocess')):
    print('Creating folder ',str(expmt_path+'/postprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/ilastik_prob').glob("mt*.tiff"): ##TO UPDATE
    path_output = Path(newpath+path_in[:-18]+'label.tif')
    list_in.append(path_in)
    list_cell.append(path_cell)
    list_out.append(path_output)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})

batch_apply(postproc_vacuole,args)

# %%
