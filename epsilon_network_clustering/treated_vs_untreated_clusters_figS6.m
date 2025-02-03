%% Load cluster (treated and untreated) and non-responder data

load("../Codon_cluster_norm_nonres_correction.mat")
load("../non_responders_per_drug_regime.mat")

%% Load clustering results
load("optimal_clusters_all_methods_CPD_scaled_102225.mat")
epsilon = methods_struct.codons.opt_epsilon;

%%
addpath('pwmetric/')
t_v_ut_pw_dist = slmetric_pw(metric_cluster_norm', untreated_metric_norm', 'eucdist');

%clustering according to the optimal epsilon value determined with the 5
%drug regime examples
cluster_array_all = opt_recluster(t_v_ut_pw_dist, epsilon, non_responders_array);
save("metric_cluster_array_perstim_ut_v_t_052324.mat", "cluster_array_all", "non_responders_array");

%% Visualize treated vs. untreated clusters (make_cluster_squares_ut_v_t.m)
num_row = 10;
num_col = 20;
subplot_order = [1:1:400];
figure_dir = '../FigS6/';
fig_path = strcat(figure_dir, 'all_drug_regime_ut_v_t_clusters');
make_cluster_squares_treated_vs_untreated(cluster_array_all, non_responders_array, num_row, num_col, subplot_order, fig_path);

%%
function cluster_array_all = opt_recluster(pw_dist, opt_epsilon, non_responders_array)

i_start = 1;
i_end = 15;
num_conditions = 15;
epsilon = opt_epsilon;
drug_regime = 1;
while i_end <= size(pw_dist, 1)
    
    pw_dist_per_drug_regime = pw_dist(i_start:i_end, :);
    non_responders_per_drug_regime = non_responders_array{drug_regime};
    cluster_array = epsilon_cluster_ut_v_t(pw_dist_per_drug_regime, epsilon, num_conditions, non_responders_per_drug_regime);
    
    cluster_array_all{drug_regime} = cluster_array;
    
    i_start = i_start + 15;
    i_end = i_end + 15;
    drug_regime = drug_regime + 1;
    
end
end