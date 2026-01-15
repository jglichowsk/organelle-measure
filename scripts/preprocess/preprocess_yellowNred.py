# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import open_organelles,neighbor_mean,batch_apply

# clean the yellow and red channels
def preprocess_yellowNred(path_in:str,path_cell:str,path_out:str,organelle:str):
    img_cell  = io.imread(str(path_cell))
    img_raw   = open_organelles[organelle](str(path_in))
    img_bkgd  = neighbor_mean(img_raw,img_cell)
    img_clean = img_raw - img_bkgd
    img_clean[img_clean<0] = 0
    img_gaussian = filters.gaussian(img_clean,sigma=0.75,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

# %%
list_in   = [] #spectral/confocal images
list_cell = [] #segmented, labelled BF images
list_out  = [] #output destination
list_orga = [] #organelle label to diff. between the three (Golgi, LD, mito)

if not os.path.exists(newpath:=Path(expmt_path+'/preprocess')):
    print('Creating folder ',str(expmt_path+'/preprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/cell_masks')).glob("GFP*.nd2"): ###TO UPDATE
    path_parts=path_in.stem.partition("_")
    date=path_parts[0]
    path_end="_".join(path_parts[1:])

    path_gl = Path(expmt_path+'/preprocess'+date)/f'gl_{path_end}.tif'
    path_mt = Path(expmt_path+'/preprocess'+date)/f'mt_{path_end}.tif' 
    path_ld = Path(expmt_path+'/preprocess'+date)/f'ld_{path_end}.tif' 

### TO UPDATE BELOW
    path_yellow = Path(str(imgs+'/'+exp+'/'+folder+'/raw'))/f"YFP_{path_cell.stem.partition('-')[2]}.nd2" 
    path_golgi = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'golgi_{path_yellow.stem.partition("_")[2]}.tif'
    path_red = Path(str(imgs+'/'+exp+'/'+folder+'/raw'))/f"red_{path_cell.stem.partition('-')[2]}_unmixed.nd2"
    path_mitochondria = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'mito_{path_red.stem.partition("_")[2][:-8]}.tif'
    path_LD = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'ld_{path_red.stem.partition("_")[2][:-8]}.tif'

    list_in.append(path_yellow)
    list_cell.append(path_cell)
    list_out.append(path_gl)
    list_orga.append("gl")

    list_in.append(path_red)
    list_cell.append(path_cell)
    list_out.append(path_mitochondria)
    list_orga.append("mt")

    list_in.append(path_red)
    list_cell.append(path_cell)
    list_out.append(path_LD)
    list_orga.append("ld")   

args = pd.DataFrame({
    "path_in": list_in,
    "path_cell": list_cell,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_yellowNred,args)