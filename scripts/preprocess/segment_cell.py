# %%
import os
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import segmentation,measure,io,util
import sys 
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from yeaz import yeaz_preprocesses,yeaz_label
from tools import load_nd2_plane,batch_apply

def segment_cells(path_in: str,path_out: str):
    img_i = load_nd2_plane(path_in,frame='yx',axes='t',idx=0)
    for prep in yeaz_preprocesses: #applies the affine transform and pixel scaling to align wide-field and confocal images. 
        img_i = prep(img_i)
    img_b = yeaz_label(img_i,min_dist=5) #applies YeaZ without GUI, using user-set min pixel distance — currently 5 px. 
    img_b = segmentation.clear_border(img_b) #clear objects touching the image border. 
    properties = measure.regionprops(img_b) #measures properties of each labelled region (here, cells) such as area, centroid, etc. 
    for prop in properties: #Omit those cells below a certain area threshold by setting them to be background.
        if prop.area < 50: # hard coded threshold, bad
            img_b[img_b==prop.label] = 0
    img_b = measure.label(img_b)
    img_o = np.zeros((512,512),dtype=int) # hard coded size, bad
    shape0,shape1 = img_b.shape
    img_o[:shape0,:shape1] = img_b

    io.imsave(str(path_out),util.img_as_uint(img_o))
    print(f"...{path_out}")
    return None

# %%
# list_in = []
# list_out = []

# imgs= master_path #path to experiment images folder
# exp= experiment_path #path to desired experiment and images
# folders= folders_list #list of experiment folders to operate on. 

# for folder in folders:
#     if not os.path.exists(newpath:=Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment'))):
#         print('Creating',str(folder+'/cell_segment'))
#         os.makedirs(newpath)
#     else:
#         print(str(folder+'/cell_segment'),'already there.')

#     for file_cell in Path(str(imgs+'/'+exp+'/'+folder+'/raw')).glob("BF*_2.nd2"): #taking the "after" BF (i.e. the one captured after spectral imaging)
#         list_in.append(file_cell)
#         file_segm = Path(newpath)/f"binCell-{file_cell.stem[3:-2]}.tif"
#         list_out.append(file_segm)
# args = pd.DataFrame({
#     "path_in":  list_in,
#     "path_out": list_out
# })

# batch_apply(segment_cells,args)
