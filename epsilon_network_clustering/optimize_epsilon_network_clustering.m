%% Objective clusters for the random drug regime examples (example_ind)

ground_truth_array{1} = {[1],[2],[3],[4],[5,6,7,8,9],[10],[11,12], [13,14,15]};
ground_truth_array{2} = {[1],[2,5,9,14,15], [3,6], [4,7],[8,13] [10], [11,12]}; % 5/9 2 14/15  
ground_truth_array{3}= {[1,2,3,13,14,15],[4,7,10],[5,8],[6],[9],[11,12]};
ground_truth_array{4} = {[1], [2], [3], [4], [5,8], [6], [7], [9,14,15], [10], [11,12], [13]}; %7 might be non responder, 14 has a big trough... mayb eby itself
ground_truth_array{5} = {[1:5,7,8,10:14], [6], [9,15]};

example_ind = [48, 91, 106, 175, 194];

%% CODONS

load("../Codon_cluster_norm_nonres_correction.mat") %metric_cluster_norm
load("../Codon_heatmap_norm_nonres_correction.mat") %metric_heatmap

%Calculate Euclidean distance between all codon feature vectors
addpath('pwmetric/')
pw_dist_nr = slmetric_pw(metric_cluster_norm', metric_cluster_norm', 'eucdist');
codons.pw_dist = pw_dist_nr;

codons.epsilon_all = [0:0.01:0.6];
codons.mis_cluster_matrix = opt_epsilon_per_method(codons.epsilon_all, pw_dist_nr, example_ind, ground_truth_array);
[codons.avg_mis_cluster, codons.min_avg_mis_cluster, codons.epsilon_choice] = get_opt_epsilon(codons.mis_cluster_matrix, codons.epsilon_all);

codons.heatmap = metric_heatmap;
%% 7 COMPONENT CPD

load("CPD_7component_cluster_scaled_0.87VE.mat") %cp_cluster_7c
load("CPD_7component_heatmap_scaled_0.87VE.mat") %cp_heatmap_7c

addpath('pwmetric/')
comp7.epsilon_all = [0:0.005:2];
pw_dist_7 = slmetric_pw(cp_cluster_7c', cp_cluster_7c', 'eucdist');
comp7.pw_dist = pw_dist_7;

comp7.mis_cluster_matrix = opt_epsilon_per_method(comp7.epsilon_all, pw_dist_7, example_ind, ground_truth_array);
[comp7.avg_mis_cluster, comp7.min_avg_mis_cluster, comp7.epsilon_choice] = get_opt_epsilon(comp7.mis_cluster_matrix, comp7.epsilon_all);

comp7.heatmap = cp_heatmap_7c;

%% 40 COMPONENT CPD

load("CPD_40component_cluster_scaled_0.95VE.mat") %cp_cluster_40c
load("CPD_40component_heatmap_scaled_0.95VE.mat") %cp_heatmap_40c

addpath('pwmetric/')
comp40.epsilon_all = [0:0.005:2.5];
pw_dist_40 = slmetric_pw(cp_cluster_40c', cp_cluster_40c', 'eucdist');
comp40.pw_dist = pw_dist_40;

comp40.mis_cluster_matrix = opt_epsilon_per_method(comp40.epsilon_all, pw_dist_40, example_ind, ground_truth_array);
[comp40.avg_mis_cluster, comp40.min_avg_mis_cluster, comp40.epsilon_choice] = get_opt_epsilon(comp40.mis_cluster_matrix, comp40.epsilon_all);

comp40.heatmap = cp_heatmap_40c;
%% 5 COMPONENT FPCA

load("fPCA_5component_cluster.mat") %fpca_cluster
load("fPCA_5component_heatmap.mat") %fpa_heatmap

addpath('pwmetric/')
pw_dist_fpca = slmetric_pw(fpca_cluster', fpca_cluster', 'eucdist');
fpca.pw_dist = pw_dist_fpca;

fpca.epsilon_all = [0:0.0001:0.08];
fpca.mis_cluster_matrix = opt_epsilon_per_method(fpca.epsilon_all, pw_dist_fpca, example_ind, ground_truth_array);
[fpca.avg_mis_cluster, fpca.min_avg_mis_cluster, fpca.epsilon_choice] = get_opt_epsilon(fpca.mis_cluster_matrix, fpca.epsilon_all);

fpca.heatmap = fpca_heatmap;

%% TRAJECTORIES
load("Traj_cluster.mat")

addpath('pwmetric/')
traj.epsilon_all = [0:0.01:0.60];
pw_dist_traj = slmetric_pw(traj_cluster', traj_cluster', 'eucdist');
traj.pw_dist = pw_dist_traj;

traj.mis_cluster_matrix = opt_epsilon_per_method(traj.epsilon_all, pw_dist_traj, example_ind, ground_truth_array);
[traj.avg_mis_cluster, traj.min_avg_mis_cluster, traj.epsilon_choice] = get_opt_epsilon(traj.mis_cluster_matrix, traj.epsilon_all);

%% Clustering across all feature spaces according to their optimal epsilon

%Choosing min optimal epsilon for now...
codons.opt_epsilon = min(codons.epsilon_choice);
codons.cluster_array_all = opt_recluster(codons.pw_dist, codons.opt_epsilon);

comp40.opt_epsilon = min(comp40.epsilon_choice);
comp40.cluster_array_all = opt_recluster(comp40.pw_dist, comp40.opt_epsilon);

comp7.opt_epsilon = min(comp7.epsilon_choice);
comp7.cluster_array_all = opt_recluster(comp7.pw_dist, comp7.opt_epsilon);

fpca.opt_epsilon = min(fpca.epsilon_choice);
fpca.cluster_array_all = opt_recluster(fpca.pw_dist, fpca.opt_epsilon);

traj.opt_epsilon = min(traj.epsilon_choice);
traj.cluster_array_all = opt_recluster(traj.pw_dist, traj.opt_epsilon);

methods_struct.codons = codons;
methods_struct.comp40 = comp40;
methods_struct.comp7 = comp7;
methods_struct.fpca = fpca; 
methods_struct.traj = traj;

save("optimal_clusters_all_methods_CPD_scaled_102225.mat", "methods_struct")
%% Run epsilon network clustering (epsilon_cluster.m) and calculate misclustering rate (min_wrong_set_match.m)
%for all 200 drug regimes and full epsilon range
function mis_cluster_matrix = opt_epsilon_per_method(epsilon_all, pw_dist, example_ind, ground_truth_array)

num_conditions = 15;
for i = 1:length(epsilon_all)
    
    i_start = 1;
    i_end = 15;
    drug_regime = 1;
    while i_end <= size(pw_dist, 2)
    
        pw_dist_per_drug_regime = pw_dist(i_start:i_end, i_start:i_end);
        cluster_array = epsilon_cluster(pw_dist_per_drug_regime, epsilon_all(i), num_conditions);
    
        cluster_array_all2{drug_regime} = cluster_array;
        
        i_start = i_start + 15;
        i_end = i_end + 15;
        drug_regime = drug_regime + 1;
    
    end
    
    for ii = 1:length(example_ind)
        
        method_cluster.clusters = cluster_array_all2{example_ind(ii)};
        method_cluster.set_index = 1:length(method_cluster.clusters);
        
        ground_truth.clusters = ground_truth_array{ii};
        ground_truth.set_index = 1:length(ground_truth.clusters);
        
        [wrong_num_min, set_matches_min] = min_wrong_set_match(ground_truth, method_cluster);
        mis_cluster_matrix(ii, i) = wrong_num_min;
    end
    

end
end

%% Find the epsilon that yields the minimum average misclustering rate
function [avg_mis_cluster, min_avg_mis_cluster, epsilon_choice] = get_opt_epsilon(mis_cluster_matrix,epsilon_all)
    avg_mis_cluster = mean(mis_cluster_matrix);
    [sorted_avg_mis_cluster, sorted_ind] = sort(avg_mis_cluster);
    [min_avg_mis_cluster, min_ind] = min(avg_mis_cluster);
    potential_epsilon = find(sorted_avg_mis_cluster == min_avg_mis_cluster);
    
    epsilon_choice = epsilon_all(sorted_ind(potential_epsilon));
end

%% Cluster each drug regime with the optimal epsilon value
function cluster_array_all = opt_recluster(pw_dist, opt_epsilon)

i_start = 1;
i_end = 15;
num_conditions = 15;
epsilon = opt_epsilon;
drug_regime = 1;
while i_end <= size(pw_dist, 2)
    
    pw_dist_per_drug_regime = pw_dist(i_start:i_end, i_start:i_end);
    cluster_array = epsilon_cluster(pw_dist_per_drug_regime, epsilon, num_conditions);
    
    cluster_array_all{drug_regime} = cluster_array;
    
    i_start = i_start + 15;
    i_end = i_end + 15;
    drug_regime = drug_regime + 1;
    
end
end