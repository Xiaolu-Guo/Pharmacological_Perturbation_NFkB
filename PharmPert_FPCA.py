import pickle, os, skfda
import numpy as np
import matplotlib.pyplot as plt
import scipy.io, scipy.stats
from skfda.exploratory.visualization import FPCAPlot
from skfda.preprocessing.dim_reduction import FPCA


data_mat = scipy.io.loadmat('../data/27-Sep-2022_NFkBsensdata.mat')

data_array = np.array(data_mat['alldata_traj_CP'])
data_array = data_array[:, :, :, 1:21, :] #Excluding untreated (index 1)

#Reorganizing into (3000, 481) matrix
cc = 0
data_matrix = np.empty((3000, data_array.shape[4]))
while cc < 3000:
    for drug_ind in np.arange(0, data_array.shape[2]):
        for drug_dose_ind in np.arange(0, data_array.shape[3]):
            for stim_ind in np.arange(0, data_array.shape[0]):
                for stim_dose_ind in np.arange(0, data_array.shape[1]):
                    traj = data_array[stim_ind, stim_dose_ind, drug_ind, drug_dose_ind, :]
                    data_matrix[cc, :] = traj

                    cc += 1

#Discretizing data matrix
data_grid = skfda.representation.grid.FDataGrid(data_matrix)

#Determining optimal number of components, using 10 component decomposition
nc = 10
fpca_discretized = FPCA(n_components = nc, centering = True)
ve = fpca_discretized.explained_variance_ratio_
ve_accum = np.empty(ve.shape)
ve_sum = 0
for i in np.arange(0, len(ve)):
    ve_sum += ve[i]
    ve_accum[i] = ve_sum

#Plotting cumulative variance explained for all 10 components
plt.plot(np.arange(0, 20), ve_accum)
plt.legend()
plt.show()

#Final decomposition using 5 components
nc = 5
fpca_discretized = FPCA(n_components = nc, centering = True)
fpca_scores = fpca_discretized.fit_transform(data_grid)
fpca_components = fpca_discretized.components_

#Saving results
results_path = '../results'
pickle.dump(fpca_scores, open(os.path.join(results_path, 'fpca_5c_scores_090623.pkl'), 'wb'))
pickle.dump(fpca_components, open(os.path.join(results_path, 'fpca_5c_components_090623.pkl'), 'wb'))