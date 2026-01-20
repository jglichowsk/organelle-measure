# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import skeletonize_zbyz,batch_apply

def postproc_ER(path_in: str,path_out: str, threshold=0.5):
    bkgdprob=io.imread(path_in)
    orgprob=1-bkgdprob
    img_in=(orgprob>threshold)
    img_ske=skeletonize_zbyz(img_in)
    io.imsave(
        str(path_out),
        util.img_as_uint(img_ske)
    )
    return None
# %%
list_in   = []
list_out  = []
list_orga = []

if not os.path.exists(newpath:=Path(expmt_path+'/postprocess')):
    print('Creating folder ',str(expmt_path+'/postprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/ilastik_prob').glob("er*.tiff"): ##TO UPDATE
    path_output=Path(newpath+path_in[:-18]+'label.tif')
    list_in.append(path_in)
    list_out.append(path_output)

args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out,
})

batch_apply(postproc_ER,args)