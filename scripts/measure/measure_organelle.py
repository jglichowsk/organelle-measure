# %%
import numpy as np
import pandas as pd
from pathlib import Path
from skimage import io,measure, util
# from organelle_measure.tools import batch_apply
# from organelle_measure.pathing_vars import master_path, experiment_path, folders_list


#List of organelle abbreviations that are used in file naming convention.
organelles = [
    "px",
    "er",
    "gl",
    "mt",
    "ld",
    "vo"
]

def parse_meta_organelle(name):
    """name is the stem of the ORGANELLE label image file."""
    #assign various labels according to file naming convention.
    organelle = name.partition("_")[0]
    labels=name.split('_')
    experiment=labels[1]
    # condition=labels[3]
    field=labels[2]     
    #then save as dictionary
    return {
        "experiment": experiment,
        # "condition":  condition,
        # "hour":       3,
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

    img_orga = io.imread(str(path_in)) #read in organelle and cell mask files
    img_cell = io.imread(str(path_cell))
    if img_cell.shape[0]>1: #if the segmentation file is a time lapse...
        img_cell=img_cell[0,:,:].astype(int) #take the first frame.

    dfs = [] #initialize list for dataframe output
    for cell in measure.regionprops(img_cell): #for each cell in the cell mask image...
        meta["idx-cell"] = cell.label #read out cell-ID
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
    print(f">>> finished {path_out}.") #send completion message. 
    return None

# %%
pi=r'C:\Users\jglic\Downloads\12-16-2025 erg6-sec61 haploid\ld_erg6-sec61_fov1_mask.tiff'
pc=r'C:\Users\jglic\Downloads\12-16-2025 erg6-sec61 haploid\20251216_BF-timelapse_cell-segm.tif'
po=r'C:\Users\jglic\Downloads\12-16-2025 erg6-sec61 haploid\20251216_BF-timelapse_org-stats.csv'
measure1organelle(pi,pc,po)
# %%
# list_in=[]; list_cell=[]; list_out=[]; #path lists for function batch process. In=organelle, cell=BF segmented tif, out=export path.

# imgs= master_path #path to experiment images folder
# exp= experiment_path #path to desired experiment and images
# folders= folders_list #list of experiment folders to operate on. 

# for folder in folders:
#     if not os.path.exists(newpath:=Path(imgs+'/'+exp+'/'+folder+'/org_measure')):
#         print('Creating folder.')
#         os.makedirs(newpath)
#     else:
#         print(str(folder+'/org_measure'),'already there.')
#     for path_cell in Path(str(imgs+'/'+exp+'/'+folder+'/cell_segment')).glob("*.tif"):
#             for organelle in organelles:
#                 path_in = Path(imgs+'/'+exp+'/'+folder+'/postprocess')/f"label-{organelle}_{path_cell.stem.partition('-')[2]}.tiff"
#                 path_out = Path(str(newpath))/f"{organelle}_{path_cell.stem.partition('-')[2]}.csv"

#                 list_in.append(path_in)
#                 list_cell.append(path_cell)
#                 list_out.append(path_out)
        
# args = pd.DataFrame({
#     "path_in":   list_in,
#     "path_cell": list_cell,
#     "path_out":  list_out
# })

# batch_apply(measure1organelle,args)

# %%
