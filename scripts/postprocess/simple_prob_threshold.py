# %%
#Apply pixel probability threshold to return binary mask of organelle signal pixels. 
import numpy as np
from pathlib import Path
from skimage import io, util

def apply_thresh(img_path:str, threshold=0.5):
    org_prob=io.imread(img_path) #read in pixel probability image
    bkgd_prob=1-org_prob
    org_mask=(org_prob>threshold) #intensity counted as true signal when greater than threshold.
    io.imsave(str(img_path)[:-18]+'label.tif', util.img_as_uint(org_mask))

# %%
# image_dir=r'C:\Users\jglic\Downloads\12-16-2025 erg6-sec61 haploid'
# for image_path in (Path(image_dir)).glob("*Probabilities.tiff"):
#     # print(image_path)
#     apply_thresh(image_path)

# %%
