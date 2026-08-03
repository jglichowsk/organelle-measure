# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,segmentation
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import batch_apply

def postproc_globular(path_in: str,path_ref: str,path_out: str, threshold=0.6):
    bkgdprob=io.imread(str(path_in))
    img_ref = io.imread(str(path_ref))
    orgprob=1-bkgdprob
    img_in=(orgprob>threshold)
    img_out=segmentation.watershed(-img_ref,mask=img_in)
    io.imsave(
        str(path_out),
        # util.img_as_ubyte(img_out)
        util.img_as_uint(img_out)
    )
    return None
# %%
list_in=[]
list_ref=[]
list_out=[]

# organelles = ["px","ld","gl"]
organelles=["ld"]

if not os.path.exists(newpath:=Path(expmt_path+'/postprocess')):
    print('Creating folder ',str(expmt_path+'/postprocess'))
    os.makedirs(newpath)

for organelle in organelles:
    for path_in in Path(expmt_path+'/ilastik_prob').glob(f"{organelle}*.tiff"):
        path_parts=path_in.stem.split("_")
        path_end="_".join(path_parts[:-1])

        path_out=Path(newpath)/f'{path_end}_label.tif'
        path_ref=Path(expmt_path+'/preprocess')/f"{path_end}.tif"
        # path_output=Path(newpath+path_in.stem[:-18]+'label.tif')
        # path_ref = (Path(expmt_path+'/preprocess'))/f"{path_in.stem[:-19]}.tif" 
        list_in.append(path_in)
        list_ref.append(path_ref)
        list_out.append(path_out)
args = pd.DataFrame({
    "path_in":  list_in,
    "path_ref": list_ref,
    "path_out": list_out
})
# %%
batch_apply(postproc_globular,args)

# %%
