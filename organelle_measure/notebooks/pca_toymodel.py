import numpy as np
from sklearn.decomposition import PCA

import torch

data = np.empty((100,2))
data[:,0] = np.arange(100)
data[:,1] = 5*np.arange(100)+12

pca = PCA(n_components=2)
pca.fit(data)

pca.components_
pca.explained_variance_ratio_
data
np.dot(data,pca.components_[1])

