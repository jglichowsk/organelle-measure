# %%
#all we're doing is rotating and rescaling the camera image to align with those from the confocal detector(s)
import numpy as np
from skimage import transform,util,io
# from tools import load_nd2_plane
# %%
#Read in camera image to be transformed.
# path_in = input() #path of image file to transform.
# path_in=r'C:\Users\jglic\Downloads\2-27-26 eyrbow\BF-frame0_sum1-7_02272026_eyrbow_glucose-2.0_fov1.tif'
# path_in=r'C:\Users\jglic\Downloads\11-7-25 tom70-tdtomato\BF_11072025_ey2795-tom70_glucose-2.0_fov1.tif'
# path_in=r'C:\Users\jglic\Downloads\2-25-26 3c-hap-selected\BF-timelapse_sum1-4_02252026_3c-hap_fov1.tif'
path_in=r'C:\Users\jglic\Downloads\SUM_BF-timelapse_03042026_eyrbow_glucose-2.0_fov2.tif'

#Rescale full size camera image down to confocal sizing if necessary.
img_2048=io.imread(str(path_in))

# params=[.9905,  .9912,  19.52, -28.11, 87.82, 6.619e-03, -5.874e-03] #1024 -> 1024 params

# params=[ 1.060,  1.060, -9.804, -25.15,  87.70, 6.480e-03, -5.700e-03] ###
# params=[ 1.055e+00,  1.055e+00, -1.000e+01, -2.400e+01,  8.770e+01, 3.492e-03, -3.600e-03]
# params=[ 1.055e+00,  1.055e+00, -1.100e+01, -2.400e+01,  8.750e+01, 0.000e+00, -1.372e-03]
params=[ 1.056e+00,  1.055e+00, -9.700e+00, -2.550e+01,  8.770e+01, 6.000e-03, -6.451e-03] #powell
scale_factor=0.25 #initial scaling; 0.5 for 1024, 0.25 for 2048


sx,sy,tx,ty,shx,shy=params[0],params[1],params[2],params[3],params[5],params[6] #unpack affine transform parameters.
theta=params[4] #unpack transform.rotate()'s parameter.
tform=transform.AffineTransform(scale=(sx,sy),translation=(tx,ty),shear=(shx,shy)) #generate the transform with given params.

if len(img_2048.shape)==2: #for time-snapshot (2D) camera images.
    img_downsc=transform.rescale(img_2048,scale_factor, anti_aliasing=True) #Downscaling dimensions of camera image.
    shape0,shape1 = img_downsc.shape #Assign downscaled image values to appropriate spots in hard coded matrix.
    img_out = np.zeros((512,512),dtype=float) #Hard coded size, bad.
    img_out[:shape0,:shape1] = img_downsc/np.max(img_downsc) #account for dimension cutoff in full camera fov. Also normalize image here.
    camera_rot=transform.rotate(img_out,theta) #rotate the raw image.
    camera_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.

elif len(img_2048.shape)==3: #for time series camera captures.
    img_downsc=transform.rescale(img_2048,scale_factor, anti_aliasing=True,channel_axis=0) #Downscaling dimensions of camera image.
    shape0,shape1,shape2 = img_downsc.shape #t,y,x; 
    img_out = np.zeros((shape0,512,512),dtype=float)
    img_out[:shape0,:shape1,:shape2] = img_downsc/np.max(img_downsc) #Assign downscaled image values to appropriate spots in hard coded matrix. Also normalize image here.
    camera_warped=[] #initalize empty list to append processed slices to.
    for i in range(img_out.shape[0]): #for each time point in series
        img_dummy=img_out[i,:,:] #create dummy array to perform operations and maintain changes. 
        camera_rot=transform.rotate(img_dummy,theta) #rotate the raw image.
        dummy_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.
        camera_warped.append(dummy_warped)
    camera_warped=np.asarray(camera_warped) #turn nested list into 3D array.
else:
    print('Outside scope of current code. Please update.') 

io.imsave(path_in[0:-4]+'_afftransf.tif',util.img_as_float32(camera_warped)) #save edited image as new file with updated name.

# %%
