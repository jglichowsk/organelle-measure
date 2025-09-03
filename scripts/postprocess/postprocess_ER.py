# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util
from organelle_measure.tools import skeletonize_zbyz,batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

### Obsolete fnc to read in h5py file type
# import h5py
# def postprocess_ER(path_in,path_out):
#     with h5py.File(str(path_in)) as f_in:
#         img_in = f_in["exported_data"][:]
#     img_in = (img_in[1]>img_in[0])
#     img_ske = skeletonize_zbyz(img_in)
#     io.imsave(
#         str(path_out),
#         util.img_as_ubyte(img_ske)
#     )
#     return None

def postproc_ER(path_in,path_out):
    bkgdprob=io.imread(path_in)
    orgprob=1-bkgdprob
    img_in=(orgprob>bkgdprob)
    img_ske=skeletonize_zbyz(img_in)
    io.imsave(
        str(path_out),
        util.img_as_ubyte(img_ske)
    )
    return None
# %%
list_in=[]; list_out=[]; #path lists for function batch process

imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 

for folder in folders:
    if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/postprocess')):
        print('Creating folder.')
        os.makedirs(newpath)
    else:
        print(str(folder+'/postprocess'),'already there.')
    for path_binary in (Path(str(imgs+'/'+exp+'/'+folder+'/ilastik prob'))).glob("ER*Probabilities.tiff"):
        path_output = Path(str(newpath))/f"label-{path_binary.stem[:-14]}.tiff"
        list_in.append(path_binary)
        list_out.append(path_output)
args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out,
})

batch_apply(postproc_ER,args)

# %%
# list_i = []
# list_o = []
# for path_binary in Path("images/preprocessed/paperRebuttal").glob("probability_ER*.h5"):
#     path_output = Path("images/labelled/paperRebuttal")/f"label-{path_binary.stem.partition('_')[2]}.tiff"
#     list_i.append(path_binary)
#     list_o.append(path_output)
# args = pd.DataFrame({
#     "path_in":  list_i,
#     "path_out": list_o
# })
# batch_apply(postprocess_ER,args)
