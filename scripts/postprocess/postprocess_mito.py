# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,util,measure
from organelle_measure.tools import batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list


def postproc_mito(path_in,path_out):
    bkgdprob=io.imread(str(path_in))
    orgprob=1-bkgdprob
    img_in=(orgprob>bkgdprob)
    img_out=measure.label(img_in)
    io.imsave(
        str(path_out),
        # util.img_as_ubyte(img_out)
        util.img_as_uint(img_out)
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
    for path_binary in (Path(str(imgs+'/'+exp+'/'+folder+'/ilastik prob'))).glob("mito*Probabilities.tiff"):
        path_output = Path(str(newpath))/f"label-{path_binary.stem[:-14]}.tiff"
        list_in.append(path_binary)
        list_out.append(path_output)
args = pd.DataFrame({
    "path_in":  list_in,
    "path_out": list_out
})

batch_apply(postproc_mito,args)
# %%
