def normalize_per_drug(tensor):
    """
        Normalizes the tensor along the drug dimension (axis = 3)

        :param tensor: tensor

        :return tensor: tensor
            Max-normalized tensor
    """

    drugs = [
    "TSA",
    "BTRCPi",
    "CHL",
    "CHX",
    "IKKi",
    "MG132",
    "PDTC",
    "PP2Ai",
    "Sel",
    "TAK1i"]

    for i_drug, _ in enumerate(drugs):
        tensor[:, :, i_drug, :, :] /= np.nanmax(tensor[:, :, i_drug, :, :])
    
    return tensor

def calculate_var_explained(tensor, rank_range):

    """
        Performs CP decomposition and calculated its error/variance explained
        for a range of ranks
        
        :param tensor: tensor
            Original, normalized, and processed tensor
        
        :param rank_range: array
            Ranks to be used in CP decomposition 
    
    """
    tensor_mask = np.isfinite(tensor)

    ve_dict = dict()
    for r in rank_range:
        sig_tensor_fact =  non_negative_parafac(tensor, init='random', rank=r, mask=tensor_mask, n_iter_max=9000, tol=1e-6, random_state=1)
        sig_tensor_fact = cp_normalize(sig_tensor_fact)

        tErr = np.nanvar(tl.cp_to_tensor(sig_tensor_fact) - tensor)
        R2X = 1.0 - (tErr / np.nanvar(tensor))

        ve_dict['Rank ' + str(r)] = R2X
    
    #Saving results
    filename = 'CPD_VE_rank_' + str(rank_range[0]) + '-' + str(rank_range[1]) + '.pkl'
    outpath = os.path.join(results_path, filename)
    pickle.dump(ve_dict, open(outpath, 'wb'))


def main():

    global results_path
    results_path = '/results'


    #Loading data matrix, converting it to a tensor
    data_mat = scipy.io.loadmat('../data/27-Sep-2022_NFkBsensdata.mat')
    data_array = np.array(data_mat['alldata_met_CP'])
    tensor = tl.tensor(data_array[:, :, :, 1:21, :]) #not including untreated condition

    #Splitting, flipping, and concatenating tensor to test decomposition's dependence on time point order
    tensor_reverse_1 = tensor[:, :, :, :, 241:482]
    tensor_reverse_2 = tensor[:, :, :, :, 0:241]
    tensor_reverse = np.concatenate((tensor_reverse_1, tensor_reverse_2), axis = 4)

    tensor_norm = normalize_per_drug(tensor)
    #tensor_norm = normalize_per_drug(tensor_reverse)

    #Pre-processing
    tensor_mask = np.isfinite(tensor_norm)
    tensor_fin = np.nan_to_num(tensor_norm)


    #Running decomposition with rank = 7
    r = 7
    sig_tensor_fact =  non_negative_parafac(tensor_fin, init='random', rank=r, mask=tensor_mask, n_iter_max=9000, tol=1e-6, random_state=1)
    sig_tensor_fact = cp_normalize(sig_tensor_fact)

    #Calculating the difference between the original and reconstructed facotrized tensors
    tErr = np.nanvar(tl.cp_to_tensor(sig_tensor_fact) - tensor_fin)

    #Normalizing and calculating the variance explained
    R2X = 1.0 - (tErr / np.nanvar(tensor))
    print('Variance explained: ' + R2X)

    #Saving results (factors and weights across all components)
    filename = 'CPD_factorized_tensor_rank_' + str(r) + '.pkl'
    outpath = os.path.join(results_path, filename)
    pickle.dump(sig_tensor_fact, open(outpath, 'wb'))

    #For determining the optimal number of components
    #calculate_var_explained(tensor_fin, np.arange(1, 9, 1))




if __name__ == "__main__":
    import numpy as np
    import tensorly as tl
    from tensorly.decomposition import non_negative_parafac, parafac
    from tensorly.cp_tensor import cp_normalize
    import scipy.io, scipy.stats
    import pickle, os