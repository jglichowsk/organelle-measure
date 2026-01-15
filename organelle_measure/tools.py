#%%
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import ndimage as ndi
from skimage import io,morphology,measure,segmentation

# BEGIN of ND2reader wrapper:
from nd2reader import ND2Reader

def get_nd2_size(path:str):
    with ND2Reader(path) as images:
        size = images.sizes
    return size

# def load_nd2_plane(path:str,frame:str='cyx',axes:str='tz',idx:int=0,channelsum=False):
#     """read an image flexibly with ND2Reader."""
#     with ND2Reader(str(path)) as images:
#         images.bundle_axes = frame
#         images.iter_axes = axes
#         if channelsum==True:
#             img=images.get_frame(idx)
#             img_summed=img.sum(axis=0)
#         else:
#             img_summed=images[idx]
#     return img_summed.squeeze()

#Rewritten to exclude boolean channelsum argument; instead do summation in organelle fncs
def load_nd2_plane(path:str,frame:str='cyx',axes:str='tz',idx:int=0):
    """read an image flexibly with ND2Reader."""
    with ND2Reader(str(path)) as images:
        images.bundle_axes = frame
        images.iter_axes = axes
        img=images[idx]
    return img
# END of ND2reader wrapper.

#%%%
# START of organelle nd2 file opener
def open_golgi(path):
    if 'c' in get_nd2_size(path):
        return np.sum(load_nd2_plane(str(path),frame="czyx",axes='t',idx=0).astype(int),axis=0,dtype=int)
    else:
        return load_nd2_plane(str(path),frame="zyx",axes='t',idx=0).astype(int)
def open_ER(path):
    if 'c' in get_nd2_size(path):
        return np.sum(load_nd2_plane(str(path),frame="czyx",axes='t',idx=0).astype(int),axis=0,dtype=int)
    else:
        return load_nd2_plane(str(path),frame="zyx",axes='t',idx=0).astype(int)
def open_mito(path):
    if Path(path).suffix == ".nd2":
        return load_nd2_plane(str(path),frame="zyx",axes='tc',idx=1).astype(int)
    else:
        img = io.imread(str(path))
        img_raw = np.zeros([int(img.shape[0]/2),*img.shape[1:]],dtype=int)
        for z in range(img_raw.shape[0]):
            img_raw[z] = img[2*z]
        return img_raw
def open_LD(path):
    if Path(path).suffix == ".nd2":
        return load_nd2_plane(str(path),frame="zyx",axes='tc',idx=0).astype(int)
    else:
        img = io.imread(str(path))
        img_raw = np.zeros([int(img.shape[0]/2),*img.shape[1:]],dtype=int)
        for z in range(img_raw.shape[0]):
            img_raw[z] = img[2*z+1]
        return img_raw
open_organelles = {
    "peroxisome":   lambda x: load_nd2_plane(str(x),frame="zyx",axes='tc',idx=0).astype(int),
    "vacuole":      lambda x: load_nd2_plane(str(x),frame="zyx",axes='tc',idx=1).astype(int),
    "ER":           open_ER,
    "golgi":        open_golgi,
    "mitochondria": open_mito,
    "LD":           open_LD
}
# END of organelle nd2 file opener

# START of bioformats wrapper
import javabridge
import bioformats
def read_spectral_img(path):
	with ND2Reader(str(path)) as nd2_img:
		size_img = nd2_img.sizes

	sample_img  = bioformats.load_image(
					str(path), 
					c=None, z=0, t=0, series=None, index=None,
					rescale=False, wants_max_intensity=False, 
					channel_names=None
	              )
	array_img = np.empty((size_img['z'],*sample_img.shape))
	for z in range(size_img['z']):
		array_img[z] = bioformats.load_image(
							str(path), 
							c=None, z=z, t=0, series=None, index=None,
							rescale=False, wants_max_intensity=False, 
							channel_names=None
	            	   )
	return array_img
# END of bioformats wrapper

# START of vacuole postprocesser

def skeletonize_zbyz(binary3d):
    """image has values [0,1] only."""
    skeletonized = np.zeros_like(binary3d)
    for z in range(len(skeletonized)):
        skeletonized[z] = morphology.skeletonize(binary3d[z])
    return skeletonized
def watershed_zbyz(skeleton3d):
    reversed3d = ~skeleton3d
    img_dist = np.zeros_like(reversed3d,dtype=float)
    img_wtsd = np.zeros_like(reversed3d,dtype=int)
    for z in range(len(reversed3d)):
        img_dist[z] = ndi.distance_transform_edt(reversed3d[z])
        img_wtsd[z] = segmentation.watershed(-img_dist[z],mask=reversed3d[z],)
    return img_wtsd
def find_complete_rings(skeleton3d):
    core3d = np.zeros_like(skeleton3d,dtype=int)
    for z in range(core3d.shape[0]):
        skeleton3d[z] = segmentation.clear_border(skeleton3d[z])
        core3d[z] = measure.label(morphology.flood_fill(~skeleton3d[z],(0,0),False,footprint=morphology.disk(1)))
    return core3d
def intersection_over_union(bool1,bool2):
    return np.count_nonzero(bool1*bool2)/np.count_nonzero(bool1+bool2)
def find_hidden_object(expand2d,watershed2d,coord):
    """
    expand2d: 2d integer image with the objects on previous plane.
    watershed: 2d watershed image on the current plane
    label: label of the object of interest
    """
    mask_last = (expand2d==expand2d[coord])
    label = watershed2d[coord]
    mask_this = (watershed2d==label)
    iou = intersection_over_union(mask_this,mask_last)
    if iou<0.5 or (np.count_nonzero(mask_this)>1.5*np.count_nonzero(mask_last)):
        return None
    else:
        return label
def better_vacuole_img(filled3d,watershed3d):
    core3d = measure.label(filled3d)
    num_z = core3d.shape[0]
    expanded3d = np.copy(core3d)
    for prop in measure.regionprops(core3d):
        # exclude very small blobs
        if prop.area<20:
            expanded3d[core3d==prop.label] = 0
            continue
        # find the range on z axis
        z_coords = prop.coords[:,0]
        z_max,z_min = z_coords.max(),z_coords.min()
        centroid = tuple(int(coo) for coo in prop.centroid[1:])
        
        # expand towards -z direction
        ref_xpnd = expanded3d[z_min]
        for z in range(z_min-1,-1,-1):
            ref_wtsd = watershed3d[z]
            chosen = find_hidden_object(ref_xpnd,ref_wtsd,centroid)
            if chosen is None:
                break
            expanded3d[z,ref_wtsd==chosen] = prop.label
            ref_xpnd = expanded3d[z]
        # expand towards +z direction
        ref_xpnd = expanded3d[z_max] 
        for z in range(z_max+1,num_z):
            ref_wtsd = watershed3d[z]
            chosen = find_hidden_object(ref_xpnd,ref_wtsd,centroid)
            if chosen is None:
                break
            expanded3d[z,ref_wtsd==chosen] = prop.label
            ref_xpnd = expanded3d[z]
    return expanded3d
# END of vacuole postprocesser


def neighbor_mean(img_orga,img_cell):
    """
    fill the regions with cells with the mean of neighbouring background.
    """
    properties = measure.regionprops(img_cell)
    img_out = np.copy(img_orga)
    for prop in properties:
        min_row, min_col, max_row, max_col = prop.bbox
        # win_row,win_col = max_row - min_row, max_col - min_col
        win_row,win_col = 2*(max_row - min_row), 2*(max_col - min_col)
        min_row = max(min_row - win_row, 0)
        min_col = max(min_col - win_col, 0)
        max_row = min(max_row + win_row, img_cell.shape[0])
        max_col = min(max_col + win_col, img_cell.shape[1])

        bbox_bin = img_cell[min_row:max_row,min_col:max_col]
        for z in range(img_orga.shape[0]):
            bbox_orga = img_orga[z,min_row:max_row,min_col:max_col]
            samples = bbox_orga[bbox_bin==0]
            background = int(samples.mean())
            idx = np.array([[z,*coo] for coo in prop.coords])
            img_out[tuple(idx.T)] = background
    return img_out

def batch_apply(func,args:pd.DataFrame):
    results = []
    errmsgs  = []
    for entry in args.iterrows():
        try:
            func(**entry[1])
            results.append(True)
            errmsgs.append("")
        except Exception as e:
            results.append(False)
            errmsgs.append(repr(e))
    args["RESULT"] = pd.Series(results)
    args["ERR_MESSAGES"] = pd.Series(errmsgs)
    return None
