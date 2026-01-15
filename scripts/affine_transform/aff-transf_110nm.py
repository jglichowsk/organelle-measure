#all we're doing is rotating and rescaling the camera image to align with those from the confocal detector(s)
import numpy as np
from skimage import transform,util,io
# from organelle_measure.tools import load_nd2_plane

#optimized affine transform params from 1024 & 512 reference images.
params=[.9905,  .9912,  19.52, -28.11, 87.82, 6.619e-03, -5.874e-03] #1024 -> 1024 params
# params=[ 9.933e-01,  9.899e-01,  1.029e+01, -2.765e+01,  8.765e+01, 4.469e-03, -3.000e-03] #512 -> 512 params
sx,sy,tx,ty,shx,shy=params[0],params[1],params[2],params[3],params[5],params[6] #unpack affine transform parameters.
theta=params[4] #unpack transform.rotate()'s parameter.
tform=transform.AffineTransform(scale=(sx,sy),translation=(tx,ty),shear=(shx,shy)) #generate the transform with given params.

# print("Affine transform:")
# print(
#     f'Scale: ({transf_affine.scale[0]:.4f}, {transf_affine.scale[1]:.4f}), '
#     f'Translation: ({transf_affine.translation[0]:.4f}, '
#     f'{transf_affine.translation[1]:.4f}), '
#     f'Rotation: {transf_affine.rotation:.4f}')

# path_in = input() #path of image file to transform.
path_in=r'C:\Users\jglic\OneDrive - Washington University in St. Louis\Documents\School\WashU\Mukherji Lab\Experiment Images\haploids\live imaging\10-31-25 erg6-sec61\20251031_124607_221__Count00000_ChannelEmpty_Seq0000.tif'
# path_in=r'C:\Users\jglic\OneDrive - Washington University in St. Louis\Documents\School\WashU\Mukherji Lab\Experiment Images\affine-transf-coords\11-2-25 110 nm beads\TRITC_512_110nm-beads_fov2-sub.tif'

img_in=io.imread(str(path_in)) #load in bright-field time series.
if len(img_in.shape)==2: #for time-snapshot (2D) camera images.
    camera_rot=transform.rotate(img_in,theta) #rotate the raw image.
    camera_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.

elif len(img_in.shape)==3: #for time series camera captures.
    camera_warped=[] #initalize empty list to append processed slices to.
    for i in range(img_in.shape[0]): #for each time point in series
        img_dummy=img_in[i,:,:] #create dummy array to perform operations and maintain changes. 
        camera_rot=transform.rotate(img_dummy,theta) #rotate the raw image.
        dummy_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.
        camera_warped.append(dummy_warped)
    camera_warped=np.asarray(camera_warped) #turn nested list into 3D array.
else:
    print('Outside scope of current code. Please update.') 

io.imsave(path_in[0:-4]+'_afftransf.tif',util.img_as_float32(camera_warped)) #save edited image as new file with updated name.






#List of processes to apply to bright-field image to match alignment and pixel sizing of confocal captures. 
# preprocesses = [
    # lambda x: transform.rotate(x,90),
    # lambda x: transform.warp(x,transf_affine.inverse)]
 #Make use of lambda function to iteratively apply these various processes to each input image. 

