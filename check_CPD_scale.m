%% Calculates a scaling factor for each component (based on their temporal 
%% pattern's AUC) in order to generate a unique decomposition

%% 7 component CPD
filename = 'tensor_traj_10drug_norm_7c_noun_allt.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
trajectories = double(tensor_factors{1, 5});
tensor_weights_7c = double(sig_fact_tensor.weights);

x = 1:481;

for i = 1:7

    AUC = trapz(x, trajectories(:, i));
    scale_all_7c(i) = 20 / AUC;
    
    AUC_new(i) = trapz(x, trajectories(:, i) * scale_all_7c(i));

end

%% 40 component CPD

filename = 'tensor_traj_10drug_norm_0.95VE_noun_allt_new.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
trajectories = double(tensor_factors{1, 5});
tensor_weights_40c = double(sig_fact_tensor.weights);

x = 1:481;

for i = 1:40

    AUC = trapz(x, trajectories(:, i));
    scale_all_40c(i) = 20 / AUC;

end

%% Save results
data_ttl = strcat("CPD_scale_7c_40c_120423", ".mat");
save(data_ttl, 'tensor_weights_7c', 'scale_all_7c', ...
               'tensor_weights_40c', 'scale_all_40c'); 