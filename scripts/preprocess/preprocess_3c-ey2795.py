# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import util,io,filters
import os
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path, red_tag, yellow_tag, blue_tag
from tools import open_organelles,neighbor_mean,batch_apply, load_nd2_plane

def open_ynb(path: str):
    """
    Args: File path to ynb nd2 capture
    Outputs: two images; one for each channel.
    """
    return load_nd2_plane(str(path),frame="zyx",axes='c',idx=0), load_nd2_plane(str(path),frame="zyx",axes='c',idx=1)
def open_red(path:str):
    """
    Args: File path to red nd2 capture
    Outputs: single image
    """
    return load_nd2_plane(str(path),frame="zyx",axes='t',idx=0)
def preprocess_red(path_in: str,path_out: str,organelle: str):
    """
    Args: File path to raw image, where to save processed image, organelle identifier
    Outputs: None; saves processed image to desginated location.
    """
    img_raw = open_red(str(path_in))
    img_gaussian = filters.gaussian(img_raw,sigma=0.4,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))
    return None
def preprocess_ynb(path_in: str,path_out: str,organelle: str):
    """
    Args: File path to raw image, where to save processed image, organelle identifier
    Outputs: None; saves processed image to desginated location.
    """
    img_vo, img_gl = open_ynb(str(path_in))
    if organelle=='vo':
        img_gaussian=filters.gaussian(img_vo,sigma=0.4,preserve_range=True).astype(int)
    else:
        img_gaussian=filters.gaussian(img_gl,sigma=0.4,preserve_range=True).astype(int)
    io.imsave(str(path_out),util.img_as_uint(img_gaussian))

# %%
list_in   = [] #spectral/confocal images
list_out  = [] #output destination
list_orga = [] #organelle label to diff. between the three (Golgi, LD, mito)

if not os.path.exists(newpath:=Path(expmt_path+'/preprocess')):
    print('Creating folder ',str(expmt_path+'/preprocess'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/raw').glob(f'{yellow_tag}*.nd2'):
    path_parts=path_in.stem.partition("_")
    path_end="".join(path_parts[1:])

    path_gl = Path(newpath)/f'gl{path_end}.tif' 
    list_in.append(path_in)
    list_out.append(path_gl)
    list_orga.append("gl")

    path_vo = Path(newpath)/f'vo{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_vo)
    list_orga.append("vo") 
args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})    
# %%
batch_apply(preprocess_ynb,args)
# %%
list_in   = [] #spectral/confocal images
list_out  = [] #output destination
list_orga = [] #organelle label to diff. between the three (Golgi, LD, mito)

for path_in in Path(expmt_path+'/raw').glob(f'{red_tag}*.nd2'):
    path_parts=path_in.stem.partition("_")
    path_end="".join(path_parts[1:])

    path_ld = Path(newpath)/f'ld{path_end}.tif'
    list_in.append(path_in)
    list_out.append(path_ld)
    list_orga.append("ld")  

args = pd.DataFrame({
    "path_in": list_in,
    "path_out": list_out,
    "organelle": list_orga
})
# %%
batch_apply(preprocess_red,args)

# %%