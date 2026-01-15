# This notebook shows the ilastik properties of a given image.
import numpy as np
from skimage import io,exposure,filters,feature

def structure_tensor_eigen(img,sigma):
    Axx,Axy,Ayy = feature.structure_tensor(filters.gaussian(img,sigma))
    eigen = np.zeros((*img.shape,2))
    eigen[...,0],eigen[...,1] = feature.structure_tensor_eigvals(Axx,Axy,Ayy)
    return eigen 
def hessian_matrix_eigen(img,sigma):
    """need testing!"""
    Helem = feature.hessian_matrix(
                filters.gaussian(img,sigma)
            )
    raise NotImplementedError

radii = [0.7,1.0,1.6,3.5,5.0,10.0]
names = [
    "_0-Gaussian-"
    "_1-LaplaceOfGaussian-",
    "_2-GaussianGradient-",
    "_3-DifferenceOfGaussian-",
    "_4-StructureTensorEigen-",
    "_5-HessianGaussianEigen-"
]
filters = [
    lambda x,y: filters.gaussian(x,sigma=y),
    lambda x,y: filters.laplace(filters.gaussian(x,sigma=y)),
    lambda x,y: filters.sobel(filters.gaussian(x,sigma=y)),
    lambda x,y: filters.difference_of_gaussians(x,low_sigma=y),
    lambda x,y: structure_tensor_eigen(x,y),
    lambda x,y: hessian_matrix_eigen(x,y)
]

# read an image
path  = ""
image = np.nan
# #############
for radius in radii:
    for j,lmbd in enumerate(filters):
        result = lmbd(image,radius)
        io.imsave(
            path.replace(
                path.partition('/')[0],
                "../test/ilastik_properties"
            ).replace(
                ".tif",
                f"{names[j]}{str(radius).replace('.','')}.tif"
            ),
            result
        )
