# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path, red_tag, yellow_tag
from tools import open_organelles,neighbor_mean,batch_apply

# clean the yellow and red channels
def preprocess_yellowNred(path_in: str,path_out: str,organelle: str):
    img_raw   = open_organelles[organelle](str(path_in))
    # img_bkgd  = neighbor_mean(img_raw,img_cell)
    # img_clean = img_raw - img_bkgd
    # img_clean[img_clean<0] = 0
    img_gaussian = filters.gaussian(img_raw,sigma=0.75,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

# %%
list_in   = [] #spectral/confocal images
list_out  = [] #output destination
list_orga = [] #organelle label to diff. between the three (Golgi, LD, mito)

if not os.path.exists(newpath:=Path(expmt_path+'/preprocess')):
    print('Creating folder ',str(expmt_path+'/preprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/raw').glob(f'{red_tag}*.nd2'):
    path_parts=path_in.stem.partition("_")
    path_end="".join(path_parts[1:])

    path_mt = Path(newpath)/f'mt{path_end[:-8]}.tif' 
    list_in.append(path_in)
    list_out.append(path_mt)
    list_orga.append("mt")

    path_ld = Path(newpath)/f'ld{path_end[:-8]}.tif'
    list_in.append(path_in)
    list_out.append(path_ld)
    list_orga.append("ld") 

for path_in in Path(expmt_path+'/raw').glob(f'{yellow_tag}*.nd2'):
    path_parts=path_in.stem.partition("_")
    path_end="".join(path_parts[1:])

    path_gl = Path(newpath)/f'gl{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_gl)
    list_orga.append("gl")  

args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})
# %%
batch_apply(preprocess_yellowNred,args)
# %%
