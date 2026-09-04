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

#imaging intervals in minutes
con_samprate = 60
cam_samprate = 5
# %%
def parse_meta_organelle(name: str):
    """
    Args: Name is the stem of the ORGANELLE label image file.
    Outptuts: Dictionary containing experiment metadata
    """
    #Unpack experiment metadata according to file naming convention. Field is fov #, and time is pre- or post-BF capture.
    tags=name.split('_')
    if tags[0]=='deconv':
        deconv, stk, time, organelle, date, strain, condition, field, probtag=name.split('_')
    else:
        stk, time, organelle, date, strain, condition, field, probtag=name.split('_')
    # field,time=field_time.split('-')
    return {
        "organelle":  organelle,
        "date":       date,
        "strain":     strain,
        "condition":  condition,
        "field":      field,
        "time":       time[-1]
    }

def measure1organelle(path_org: str, path_cell: str, path_out: str, metadata=None):
    """
    Args: File paths for organelle mask frames, cell mask file, and output save location.
    Outputs: Saves designated metrics to csv file at the location of path_out.
    """
    # parse metadata from filename
    if metadata is None:
        meta = parse_meta_organelle(path_org.stem)
    else:
        meta = metadata

    #read in organelle and cell mask files
    img_orga = io.imread(path_org) 
    img_cell = io.imread(path_cell)

    if img_cell.shape[0]>1: #if timelapse, choose appropriate mask frame and store abs_time info
        frame=int(meta['time'])
        abs_time = (frame-1) * con_samprate
        meta['abs_time']=abs_time
        cam_frame=int(abs_time/cam_samprate)
        img_cell=img_cell[cam_frame,:,:].astype(int)       

    dfs = [] #initialize list for dataframe output
    for cell in measure.regionprops(img_cell): #for each cell in the cell mask image...
        meta["idx_cell"] = cell.label #read out cell-ID
        min_row, min_col, max_row, max_col = cell.bbox #extract bounding box details
        img_cell_crop = cell.image #cell image mask

        if len(img_orga.shape)==3: #3D
            img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col] #crop organelle image to just include given cell
            for z in range(img_orga_crop.shape[0]): #for each z-slice in organelle image stack...
                img_orga_crop[z] = img_orga_crop[z]*img_cell_crop #apply the midplane cell segmentation mask.
        else: #2D
            img_orga_crop = img_orga[min_row:max_row,min_col:max_col]
            img_orga_crop = img_orga_crop*img_cell_crop
        
        if not meta["organelle"] == "vo": #if not vacuole, read out the following properties
            measured_orga = measure.regionprops_table(
                img_orga_crop,
                properties=('label','area','bbox_area','bbox')
            )
        else: #if vacuole...
            cell_minor_axis=cell["axis_minor_length"]
            measured_orga = measure.regionprops_table(img_orga_crop, properties=('label','area','bbox_area','bbox'))
            vo_area=measured_orga["area"]
            vo_vol_est = vo_area * cell_minor_axis
            measured_orga["area"]=vo_vol_est

        result = meta | measured_orga #join the metadata and extracted organelle metrics.
        dfs.append(pd.DataFrame(result)) #append into dataframe
    
    if len(dfs) == 0: #in case of no cells...
        print(f">>> {path_out} has no cells, skipped.") #send error message.
        return None
    
    df_orga = pd.concat(dfs,ignore_index=True) #join all single-cell dataframes outputted. 
    df_orga.rename(columns={'label':'idx-orga',"area":"volume-pixel",'bbox_area':'volume-bbox'},inplace=True) #some column labelling
    # df_orga.to_csv(str(path_out),index=False) #save dataframe as csv file
    # print(f"Finished {path_out}.") #send completion message. 
    return df_orga
# %% Generate variables to run through functions above
list_in   = [] #Initialize lists for org masks, cell masks, and output paths respectively. 
list_cell = []
list_out  = []

if not os.path.exists(newpath:=Path(expmt_path+'/org_measure')):
    print('Creating folder ',str(expmt_path+'/org_measure'))
    os.makedirs(newpath)

# organelles = ["px","er","gl","mt","ld","vo"]
organelles=['ld', 'gl', 'vo']
for organelle in organelles:
    for path_in in Path(expmt_path+'/postprocess').glob("*label.tif"): 
        meta=parse_meta_organelle(path_in.stem)
        if meta['organelle']==organelle:
            path_parts=path_in.stem.split("_")
            if path_parts[0]=='deconv':
                path_end="_".join(path_parts[4:-1])
            else:
                path_end="_".join(path_parts[3:-1])
            path_cell=Path(expmt_path+'/cell_segment')/f"BF-timelapse_{path_end}_segm-afftransf-expand.tif"
            path_out=Path(newpath)/f"{organelle}_{path_end}.csv"

            list_in.append(path_in)
            list_cell.append(path_cell)
            list_out.append(path_out)

args = pd.DataFrame({
    "path_in":   list_in,
    "path_cell": list_cell,
    "path_out":  list_out
})

# %% batch apply
# orgs=['ld', 'gl']
orgs=["ld"]
# orgs=["vo"]

for organelle in orgs:
    org_dfs=[]
    for path in list_in:
        meta=parse_meta_organelle(path.stem)
        if meta['organelle']==organelle:
            ind=list_in.index(path)
            org_df=measure1organelle(list_in[ind], list_cell[ind], list_out[ind])
            org_dfs.append(org_df)
    output_df=pd.concat(org_dfs, ignore_index=True)
    output_df.to_csv(list_out[ind],index=False)

# %% previous iteration of measure1organelle
# def measure1organelle(path_org: str, path_cell: str, path_out: str, metadata=None):
#     """
#     Args: File paths for organelle mask, cell mask, and output save location.
#     Outputs: Saves designated metrics to csv file at the location of path_out.
#     """
#     # parse metadata from filename
#     name = Path(path_org).stem
#     if metadata is None:
#         meta = parse_meta_organelle(name)
#     else:
#         meta = metadata

#     #read in organelle and cell mask files
#     img_orga = io.imread(str(path_org)) 
#     img_cell = io.imread(str(path_cell))

#     if img_cell.shape[0]>1: #if the cell segmentation file is a time lapse...
#         if meta['time']=='pre':
#             img_cell=img_cell[0,:,:].astype(int) #take the first frame for overlay.
#         elif meta['time']=='post':
#             img_cell=img_cell[-1,:,:].astype(int) #take the last frame for overlay.
#         else:
#             print('File tag not recognized. Please verify that the input masks have the correct naming convention.')        

#     dfs = [] #initialize list for dataframe output
#     for cell in measure.regionprops(img_cell): #for each cell in the cell mask image...
#         meta["idx_cell"] = cell.label #read out cell-ID
#         min_row, min_col, max_row, max_col = cell.bbox #extract bounding box details
#         img_orga_crop = img_orga[:,min_row:max_row,min_col:max_col] #crop organelle image to just include given cell
#         img_cell_crop = cell.image #same crop but for cell image
#         for z in range(img_orga_crop.shape[0]): #for each z-slice in organelle image stack...
#             img_orga_crop[z] = img_orga_crop[z]*img_cell_crop #apply the midplane cell segmentation mask.
#         if not meta["organelle"] == "vo": #if not vacuole, read out the following properties
#             measured_orga = measure.regionprops_table(
#                 img_orga_crop,
#                 properties=('label','area','bbox_area','bbox')
#             )
#         else: #if vacuole...
#             vo_area = 0 #initialize the following metrics at 0
#             vo_bbox_area = 0
#             bbox0,bbox1,bbox2,bbox3,bbox4,bbox5 = 0,0,0,0,0,0
#             for z in range(img_orga_crop.shape[0]): #for each z-slice in cell image stack...
#                 vo = measure.regionprops_table(img_orga_crop[z], properties=('label','area','bbox_area','bbox'))

#                 if len(vo["area"]) == 0: #if the vacuole signal in given z-slice is zero...
#                     continue #do nothing?
#                 if (maxblob:=max(vo["area"])) > vo_area: #run through to find max vo area
#                     vo_area = maxblob #update variable
#                     idxblob = np.argmax(vo["area"]) #record index (z-slice) of max vacuole area.
#                     vo_bbox_area = vo["bbox_area"][idxblob] #extract bbox area for that slice. 
#                     bbox0,bbox3 = z,z #rewriting some bbox params according to where max area slice is located?
#                     bbox1,bbox2,bbox4,bbox5 = [vo[f"bbox-{i}"][idxblob] for i in range(4)]
#             if vo_area==0: #if still 0 area...
#                 continue #go ahead
#             #Read out the following metrics for vacuoles.
#             measured_orga = {
#                 'label': [0],
#                 'area':  [vo_area],
#                 "bbox_area": [vo_bbox_area],
#                 "bbox-0": [bbox0],
#                 "bbox-1": [bbox1],
#                 "bbox-2": [bbox2],
#                 "bbox-3": [bbox3],
#                 "bbox-4": [bbox4],
#                 "bbox-5": [bbox5],
#             }
#         result = meta | measured_orga #join the metadata and extracted organelle metrics. 
#         dfs.append(pd.DataFrame(result)) #append into dataframe
#     if len(dfs) == 0: #in case of no cells...
#         print(f">>> {path_out} has no cells, skipped.") #send error message.
#         return None
#     df_orga = pd.concat(dfs,ignore_index=True) #join all single-cell dataframes outputted. 
#     df_orga.rename(columns={'label':'idx-orga',"area":"volume-pixel",'bbox_area':'volume-bbox'},inplace=True) #some column labelling
#     df_orga.to_csv(str(path_out),index=False) #save dataframe as csv file
#     print(f"Finished {path_out}.") #send completion message. 
#     return None
# %%
