# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skimage import io, segmentation, measure
from pathlib import Path
import os
import sys
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path

# %%
def remove_interp_pixels(img: np.ndarray):
    """
    Args: 3D Cell mask image (t,y,x) with non-integer, interpolated pixels.
    Outputs: Cell mask image with said pixel values set to zero.
    """
    nonints= (img%1!=0)
    img[nonints]=0
    return img 

for img_affpath in Path(expmt_path+'/cell_segment').glob('*_segm-afftransf.tif'):
    print(f">>>> Processing {img_affpath} <<<<")
    img_aff=io.imread(img_affpath)
    img_nobord=remove_interp_pixels(img_aff)
    io.imsave(Path(img_affpath).parent/Path(Path(img_affpath).stem+'-nobord.tif'),img_nobord)

    #expand twice in succession to better divvy up bud neck region between mother and daughter
    for t in range(img_nobord.shape[0]):
        img_nobord[t,:,:]=segmentation.expand_labels(img_nobord[t,:,:], distance=1) 
        img_nobord[t,:,:]=segmentation.expand_labels(img_nobord[t,:,:], distance=1) 
        img_nobord[t,:,:]=segmentation.expand_labels(img_nobord[t,:,:], distance=1) 

    io.imsave(Path(img_affpath).parent/Path(Path(img_affpath).stem+'-expand.tif'), img_nobord)


# %%
