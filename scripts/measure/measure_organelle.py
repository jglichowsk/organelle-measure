# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,measure
from organelle_measure.tools import batch_apply
from organelle_measure.pathing_vars import master_path, experiment_path, folders_list

def parse_meta_organelle(name):
    """name is the stem of the ORGANELLE label image file."""
    #assign various labels according to file naming convention.
    organelle = name.partition("-")[2].partition("_")[0]
    labels=name.split('_')
    experiment=labels[2]
    condition=labels[3]
    field=labels[4]     
    return {
        "experiment": experiment,
        "condition":  condition,
        "hour":       3,
        "field":      field,
        "organelle":  organelle
    }

def measure1organelle(path_in,path_cell,path_out,metadata=None):
    # parse metadata from filename
    name = Path(path_in).stem
    if metadata is None:
        meta = parse_meta_organelle(name)
    else:
        meta = metadata

    img_orga = io.imread(str(path_in))
    img_cell = io.imread(str(path_cell))
    
    dfs = []
    for cell in measure.regionprops(img_cell):
        meta["idx-cell"] = cell.label
        min_row, min_col, max_row, max_col = cell.bbox
        img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col]
        img_cell_crop = cell.image
        for z in range(img_orga_crop.shape[0]):
            img_orga_crop[z] = img_orga_crop[z]*img_cell_crop
        if not meta["organelle"] == "vacuole":
            measured_orga = measure.regionprops_table(
                img_orga_crop,
                properties=('label','area','bbox_area','bbox')
            )
        else:
            vacuole_area = 0
            vacuole_bbox_area = 0
            bbox0,bbox1,bbox2,bbox3,bbox4,bbox5 = 0,0,0,0,0,0
            for z in range(img_orga_crop.shape[0]):
                vacuole = measure.regionprops_table(
                    img_orga_crop[z],
                    properties=('label','area','bbox_area','bbox')
                )
                if len(vacuole["area"]) == 0:
                    continue
                if (maxblob:=max(vacuole["area"])) > vacuole_area:
                    vacuole_area = maxblob
                    idxblob = np.argmax(vacuole["area"])
                    vacuole_bbox_area = vacuole["bbox_area"][idxblob]
                    bbox0,bbox3 = z,z
                    bbox1,bbox2,bbox4,bbox5 = [vacuole[f"bbox-{i}"][idxblob] for i in range(4)]
            if vacuole_area==0:
                continue
            measured_orga = {
                'label': [0],
                'area':  [vacuole_area],
                "bbox_area": [vacuole_bbox_area],
                "bbox-0": [bbox0],
                "bbox-1": [bbox1],
                "bbox-2": [bbox2],
                "bbox-3": [bbox3],
                "bbox-4": [bbox4],
                "bbox-5": [bbox5],
            }
        result = meta | measured_orga
        dfs.append(pd.DataFrame(result))
    if len(dfs) == 0:
        print(f">>> {path_out} has no cells, skipped.")
        return None
    df_orga = pd.concat(dfs,ignore_index=True)
    df_orga.rename(columns={'label':'idx-orga',"area":"volume-pixel",'bbox_area':'volume-bbox'},inplace=True)
    df_orga.to_csv(str(path_out),index=False)
    print(f">>> finished {path_out.stem}.")
    return None

organelles = [
    "peroxisome",
    "ER",
    "golgi",
    "mito",
    "LD",
    "vacuole"
]

# %%
list_in=[]; list_cell=[]; list_out=[]; #path lists for function batch process. In=organelle, cell=BF segmented tif, out=export path.

imgs= master_path #path to experiment images folder
exp= experiment_path #path to desired experiment and images
folders= folders_list #list of experiment folders to operate on. 

for folder in folders:
    if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/org_measure')):
        print('Creating folder.')
        os.makedirs(newpath)
    else:
        print(str(folder+'/org_measure'),'already there.')
    for path_cell in Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment')).glob("*.tif"):
            for organelle in organelles:
                path_in = Path(imgs+'/'+exp+'/'+folder+'/postprocess')/f"label-{organelle}_{path_cell.stem.partition('-')[2]}.tiff"
                path_out = Path(str(newpath))/f"{organelle}_{path_cell.stem.partition('-')[2]}.csv"

                list_in.append(path_in)
                list_cell.append(path_cell)
                list_out.append(path_out)
        
args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})

batch_apply(measure1organelle,args)

# %%
