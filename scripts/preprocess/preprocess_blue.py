# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import open_organelles,neighbor_mean,batch_apply

def preprocess_blue(path_in: str,path_out: str,organelle: str):
    """
    Args: File path to raw image, where to save processed image, organelle identifier
    Outputs: None; saves processed image to desginated location.
    """
    img_raw=open_organelles[organelle](str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.75,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None

list_in   = [] #spectral/confocal images
list_out  = [] #output destination
list_orga = [] #organelle label to differentiate between the peroxisome and vacuoles

if not os.path.exists(newpath:=Path(expmt_path+'/preprocess')):
    print('Creating folder ',str(expmt_path+'/preprocess'))
    os.makedirs(newpath)

for path_in in (Path(expmt_path+'/raw')).glob("CFP*unmixed.nd2"): ###TO UPDATE
    path_parts=path_in.stem.partition("_")
    date=path_parts[0]
    path_end="_".join(path_parts[1:])

    path_px = Path(expmt_path+'/preprocess'+date)/f'px_{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_px)
    list_orga.append("px")

    path_vo = Path(expmt_path+'/preprocess'+date)/f'vo_{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_vo)
    list_orga.append("vo")
args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})

batch_apply(preprocess_blue,args)

# %%