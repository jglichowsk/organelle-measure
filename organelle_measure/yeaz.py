# BEGIN segment cells with YeaZ
import numpy as np
from .unet import neural_network,segment
from skimage import transform,exposure,util
import torch

#matrix characerizing the affine transform needed to align the wide-field and confocal output images.
matrix_affine = np.array(
    [[ 9.48674270e-01,  4.20567174e-02,  2.11683860e-02],
     [-4.21646062e-02,  9.50886355e-01,  2.58659871e+01],
     [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
transf_affine = transform.AffineTransform(matrix=matrix_affine)

#list of porcesses to apply to bright-field image to match alignment and pixel sizing of confocal captures. 
yeaz_preprocesses = [
    lambda x: transform.rotate(x,90),
    lambda x: transform.rescale(x,0.25),
    lambda x: transform.warp(x,transf_affine),
    lambda x: (x-x.min())/(x.max()-x.min())
] #makes use of lambda function to iteratively apply these various processes to each input image. 

def yeaz_label(img_i,min_dist,**kwargs):
    """
    Use YeaZ without GUI to segment cells in an image.
    INPUT:  
        img_i, a 2-D float image containing cells to segment. 
        min_dist, the minimum distance passed to segment()
    OUTPUT: img_o, a 2-D uint label image, same shape as img_i.
    """
    img_exposure  = exposure.equalize_adapthist(img_i)
    img_predict   = neural_network.prediction(img_exposure,mic_type="bf",device=torch.device("cpu"),**kwargs)
    img_threshold = neural_network.threshold(img_predict)
    img_segment   = segment.segment(
                                    img_threshold,img_predict,
                                    min_distance=min_dist,
                                    topology=None
                                   )
    print(img_segment.max())
    img_o = util.img_as_uint(img_segment.astype(int))
    return img_o
# END segment cells with YeaZ