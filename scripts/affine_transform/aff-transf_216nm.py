#all we're doing is rotating and rescaling the camera image to align with those from the confocal detector(s)
import numpy as np
from skimage import transform,util,io
# from organelle_measure.tools import load_nd2_plane

# path_in = input() #path of image file to transform.
# path_in=r'C:\Users\jglic\Downloads\1-22-26 myo1-mLemon\BF_01222026_myo1_glucose-2.0_fov1.tif'
path_in=r'C:\Users\jglic\Downloads\1-22-26 myo1-mLemon\BF_01222026_myo1_glucose-2.0_fov1_slice-16.tif'
img_in=io.imread(str(path_in)) #load in bright-field time series.

#optimized affine transform params from 1024 & 512 reference images.
params=[1.007, 1.008, 6.174, -22.54, 87.79, 6.487e-03, -5.678e-03] #1024 -> 512 params
sx,sy,tx,ty,shx,shy=params[0],params[1],params[2],params[3],params[5],params[6] #unpack affine transform parameters.
theta=params[4] #unpack transform.rotate()'s parameter.
tform=transform.AffineTransform(scale=(sx,sy),translation=(tx,ty),shear=(shx,shy)) #generate the transform with given params.

# print("Affine transform:")
# print(
#     f'Scale: ({transf_affine.scale[0]:.4f}, {transf_affine.scale[1]:.4f}), '
#     f'Translation: ({transf_affine.translation[0]:.4f}, '
#     f'{transf_affine.translation[1]:.4f}), '
#     f'Rotation: {transf_affine.rotation:.4f}')

if len(img_in.shape)==2: #for time-snapshot (2D) camera images.
    camera_rot=transform.rotate(img_in,theta) #rotate the raw image.
    camera_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.
    camera_warped=camera_warped/np.max(camera_warped) #normalize

elif len(img_in.shape)==3: #for time series camera captures.
    camera_warped=[] #initalize empty list to append processed slices to.
    for i in range(img_in.shape[0]): #for each time point in series
        img_dummy=img_in[i,:,:] #create dummy array to perform operations and maintain changes. 
        img_downsc=transform.rescale(img_dummy,0.5, anti_aliasing=True)
        img_out = np.zeros((512,512),dtype=float) # hard coded size, okay.
        shape0,shape1 = img_downsc.shape
        img_out[:shape0,:shape1] = img_downsc #account for dimension cutoff in full camera fov. Also normalize image here.
        camera_rot=transform.rotate(img_out,theta) #rotate the raw image.
        dummy_warped=transform.warp(camera_rot,tform.inverse) #apply the affine transform to the rotated camera image.
        camera_warped.append(dummy_warped/np.max(dummy_warped))
    camera_warped=np.asarray(camera_warped) #turn nested list into 3D array.
else:
    print('Outside scope of current code. Please update.') 

io.imsave(path_in[0:-4]+'_afftransf.tif',util.img_as_float32(camera_warped)) #save edited image as new file with updated name.