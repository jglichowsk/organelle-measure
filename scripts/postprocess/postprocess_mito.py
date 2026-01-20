# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,measure
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import batch_apply

def postproc_mito(path_in: str,path_out: str, threshold=0.5):
    bkgdprob=io.imread(str(path_in))
    orgprob=1-bkgdprob
    img_in=(orgprob>threshold)
    img_out=measure.label(img_in)
    io.imsave(
        str(path_out),
        util.img_as_uint(img_out)
    )
    return None
# %%
list_in=[]
list_out=[]

if not os.path.exists(newpath:=Path(expmt_path+'/postprocess')):
    print('Creating folder ',str(expmt_path+'/postprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/ilastik_prob').glob("mt*.tiff"): ##TO UPDATE
    path_output = Path(newpath+path_in[:-18]+'label.tif')
    list_in.append(path_in)
    list_out.append(path_output)
args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out
})

batch_apply(postproc_mito,args)
# %%
