# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import segmentation,measure,io,util
from organelle_measure.yeaz import yeaz_preprocesses,yeaz_label
from organelle_measure.tools import load_nd2_plane,batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

def segment_cells(path_in,path_out):
    img_i = load_nd2_plane(str(path_in),frame='yx',axes='t',idx=0)
    for prep in yeaz_preprocesses:
        img_i = prep(img_i)
    img_b = yeaz_label(img_i,min_dist=5)
    img_b = segmentation.clear_border(img_b)
    properties = measure.regionprops(img_b)
    for prop in properties:
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
list_in = []
list_out = []

imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 

for folder in folders:
    if not os.path.exists(newpath:=Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment'))):
        print('Creating',str(folder+'/cell_segment'))
        os.makedirs(newpath)
    else:
        print(str(folder+'/cell_segment'),'already there.')

    for file_cell in Path(str(imgs+'/'+exp+'/'+folder+'/raw')).glob("BF*_2.nd2"): #taking the "after" BF (i.e. the one captured after spectral imaging)
        list_in.append(file_cell)
        file_segm = Path(newpath)/f"binCell-{file_cell.stem[3:-2]}.tif"
        list_out.append(file_segm)
args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out
})

batch_apply(segment_cells,args)

# %%

# list_in = []
# list_out = []

# #path to experiment images folder
# imgs_path=f"C:/Users/jglic/OneDrive - Washington University in St. Louis/Documents/School/WashU/Mukherji Lab/Experiment Images"
# exp_path=f"rbow knockouts/BF_only" #path to desired experiment and images
# folder = f"7-23-24"
# if not os.path.exists(newpath:=Path(str(imgs_path+'/'+exp_path+'/'+folder+'/cell_segment'))):
#     print('Creating',str(folder,'/cell_segment'))
#     os.makedirs(newpath)
# else:
#     print(str(folder+'/cell_segment'),'already there.')
# for file_cell in Path(str(imgs_path+'/'+exp_path+'/'+folder+'/raw')).glob("BF*.nd2"):
#     list_in.append(file_cell)
#     file_segm = Path(newpath)/f"binCell-{file_cell.stem[3:]}.tif"
#     list_out.append(file_segm)
# args = pd.DataFrame({
#     "path_in":  list_in,
#     "path_out": list_out
# })
# batch_apply(segment_cells,args)

# %%
