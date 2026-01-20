# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import open_organelles,neighbor_mean,batch_apply

def preprocess_green(path_in: str,path_out: str,organelle: str):
    img_raw   = open_organelles[organelle](str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.3,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

list_in   = []
list_out  = []
list_orga = []

if not os.path.exists(newpath:=Path(expmt_path+'/preprocess')):
    print('Creating folder ',str(expmt_path+'/preprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/raw').glob("GFP*.nd2"): ##TO UPDATE
    path_parts=path_in.stem.partition("_")
    path_end="_".join(path_parts[1:])
    path_er = Path(newpath)/f'er_{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_er)
    list_orga.append("er")

args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})
batch_apply(preprocess_green,args)
# %%
