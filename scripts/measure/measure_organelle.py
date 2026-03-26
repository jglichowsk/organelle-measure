# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,measure, util
import os
import sys
sys.path.append(r'C:\Users\jglic\Documents\School\WashU\Mukherji Lab\organelle-measure\organelle_measure')
from pathing_variables import expmt_path
from tools import batch_apply

#List of organelle abbreviations that are used in file naming convention.
organelles = [
    "px",
    "er",
    "gl",
    "mt",
    "ld",
    "vo"
]

def parse_meta_organelle(name: str):
    """
    Args: Name is the stem of the ORGANELLE label image file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention. Field is fov #, and time is pre- or post-BF capture.
    organelle,date,strain,condition,field_time,label=name.split('_')
    field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
        "time":       time
    }

def measure1organelle(path_org: str,path_cell: str,path_out: str,metadata=None):
    """
    Args: File paths for organelle mask, cell mask, and output save location.
    Outputs: Saves designated metrics to csv file at the location of path_out.
    """
    # parse metadata from filename
    name = Path(path_org).stem
    if metadata is None:
        meta = parse_meta_organelle(name)
    else:
        meta = metadata

    #read in organelle and cell mask files
    img_orga = io.imread(str(path_org)) 
    img_cell = io.imread(str(path_cell))

    if img_cell.shape[0]>1: #if the cell segmentation file is a time lapse...
        if meta['time']=='pre':
            img_cell=img_cell[0,:,:].astype(int) #take the first frame for overlay.
        elif meta['time']=='post':
            img_cell=img_cell[-1,:,:].astype(int) #take the last frame for overlay.
        else:
            print('File tag not recognized. Please verify that the input masks have the correct naming convention.')        

    dfs = [] #initialize list for dataframe output
    for cell in measure.regionprops(img_cell): #for each cell in the cell mask image...
        meta["idx_cell"] = cell.label #read out cell-ID
        min_row, min_col, max_row, max_col = cell.bbox #extract bounding box details
        img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col] #crop organelle image to just include given cell
        img_cell_crop = cell.image #same crop but for cell image
        for z in range(img_orga_crop.shape[0]): #for each z-slice in organelle image stack...
            img_orga_crop[z] = img_orga_crop[z]*img_cell_crop #apply the midplane cell segmentation mask.
        if not meta["organelle"] == "vo": #if not vacuole, read out the following properties
            measured_orga = measure.regionprops_table(
                img_orga_crop,
                properties=('label','area','bbox_area','bbox')
            )
        else: #if vacuole...
            vo_area = 0 #initialize the following metrics at 0
            vo_bbox_area = 0
            bbox0,bbox1,bbox2,bbox3,bbox4,bbox5 = 0,0,0,0,0,0
            for z in range(img_orga_crop.shape[0]): #for each z-slice in organelle image stack...
                #read out the following properties to table.
                vo = measure.regionprops_table(
                    img_orga_crop[z],
                    properties=('label','area','bbox_area','bbox')
                )
                if len(vo["area"]) == 0: #if the vacuole signal in given z-slice is zero...
                    continue #go with it?
                if (maxblob:=max(vo["area"])) > vo_area: #if the max area of previously-measured slices is greater than
                    #the current slice, then...
                    vo_area = maxblob #assign that max value to current area?
                    idxblob = np.argmax(vo["area"]) #record index (z-slice) of max vacuole area.
                    vo_bbox_area = vo["bbox_area"][idxblob] #extract bbox area for that slice. 
                    bbox0,bbox3 = z,z #rewriting some bbox params according to where max area slice is located?
                    bbox1,bbox2,bbox4,bbox5 = [vo[f"bbox-{i}"][idxblob] for i in range(4)]
            if vo_area==0: #if still 0 area...
                continue #go ahead
            #Read out the following metrics for vacuoles.
            measured_orga = {
                'label': [0],
                'area':  [vo_area],
                "bbox_area": [vo_bbox_area],
                "bbox-0": [bbox0],
                "bbox-1": [bbox1],
                "bbox-2": [bbox2],
                "bbox-3": [bbox3],
                "bbox-4": [bbox4],
                "bbox-5": [bbox5],
            }
        result = meta | measured_orga #join the metadata and extracted organelle metrics. 
        dfs.append(pd.DataFrame(result)) #append into dataframe
    if len(dfs) == 0: #in case of no cells...
        print(f">>> {path_out} has no cells, skipped.") #send error message.
        return None
    df_orga = pd.concat(dfs,ignore_index=True) #join all single-cell dataframes outputted. 
    df_orga.rename(columns={'label':'idx-orga',"area":"volume-pixel",'bbox_area':'volume-bbox'},inplace=True) #some column labelling
    df_orga.to_csv(str(path_out),index=False) #save dataframe as csv file
    print(f"Finished {path_out}.") #send completion message. 
    return None
# %%

list_in   = [] #Initialize lists for org masks, cell masks, and output paths respectively. 
list_cell = []
list_out  = []

if not os.path.exists(newpath:=Path(expmt_path+'/org_measure')):
    print('Creating folder ',str(expmt_path+'/org_measure'))
    os.makedirs(newpath)

for path_in in Path(expmt_path+'/postprocess').glob('*label.tif'):
    path_parts=path_in.stem.split("_")
    path_end="_".join(path_parts[:-1])
    cell_parts=path_in.stem.split('-')
    cell_end="-".join(cell_parts[:2])[3:]

    path_cell=Path(expmt_path+'/cell_segment')/f"BF-timelapse_{cell_end}_segm.tif"
    path_out=Path(newpath)/f'{path_end}.csv'

    list_in.append(path_in)
    list_cell.append(path_cell)
    list_out.append(path_out)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})

# %%
# batch_apply(measure1organelle,args)
for i in range(len(list_in)):
    measure1organelle(list_in[i],list_cell[i],list_out[i])

# %%
