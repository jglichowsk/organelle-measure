# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
from organelle_measure.tools import open_organelles,neighbor_mean,batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

# clean the yellow and red channels
def preprocess_yellowNred(path_in,path_cell,path_out,organelle):
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
list_cell = [] #segmented, labelled BF iamges
list_out  = [] #output destination
list_orga = [] #organelle label to diff. between the three (Golgi, LD, mito)

imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 
print("Creating folder as needed.")
for folder in folders:
    if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/preprocess')):
        print('Creating',str(folder+'/preprocess'))
        os.makedirs(newpath)
    else:
        print(str(folder+'/preprocess'),'already there.')
    for path_cell in (Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment'))).glob("*.tif"):
        path_yellow = Path(str(imgs+'/'+exp+'/'+folder+'/raw'))/f"YFP_{path_cell.stem.partition('-')[2]}.nd2"
        path_golgi = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'golgi_{path_yellow.stem.partition("_")[2]}.tif'
        
        list_in.append(path_yellow)
        list_cell.append(path_cell)
        list_out.append(path_golgi)
        list_orga.append("golgi")

        path_red = Path(str(imgs+'/'+exp+'/'+folder+'/raw'))/f"red_{path_cell.stem.partition('-')[2]}_unmixed.nd2"
        path_mitochondria = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'mito_{path_red.stem.partition("_")[2][:-8]}.tif'

        list_in.append(path_red)
        list_cell.append(path_cell)
        list_out.append(path_mitochondria)
        list_orga.append("mitochondria")

        path_LD = Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))/f'LD_{path_red.stem.partition("_")[2][:-8]}.tif'

        list_in.append(path_red)
        list_cell.append(path_cell)
        list_out.append(path_LD)
        list_orga.append("LD")   
args = pd.DataFrame({
    "path_in": list_in,
    "path_cell": list_cell,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_yellowNred,args)

# %%