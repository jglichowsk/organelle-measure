# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
from organelle_measure.tools import open_organelles,neighbor_mean,batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

def preprocess_green(path_in,path_out,organelle):
    img_raw   = open_organelles[organelle](str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.3,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

# %%
list_in   = []
list_out  = []
list_orga = []

imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 

print("Creating folders as needed.")
for folder in folders:
    if not os.path.exists(newpath:=Path(str(imgs+'/'+exp+'/'+folder+'/preprocess'))):
        print('Creating',str(folder+'/preprocess'))
        os.makedirs(newpath)
    else:
        print(str(folder+'/preprocess'),'already there.')

    for path_in in Path(str(imgs+'/'+exp+'/'+folder+'/raw')).glob("GFP*.nd2"):
        path_ER = Path(newpath)/f'ER_{path_in.stem.partition("_")[2]}.tif'
        list_in.append(path_in)
        list_out.append(path_ER)
        list_orga.append("ER")

args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})
batch_apply(preprocess_green,args)
# %%
