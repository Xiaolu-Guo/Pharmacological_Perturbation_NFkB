Description of the codes:

(1) pharm_pert_sens_analysis/run_pharm_pert_sens_analysis.m: run simulation to get 3000 trajectories ( NFkB dynamics response to 15 stimuli under 10 drug  across 20 doses for each drug) [uploaded]. 
	*ODE model can be accessed at ___

(2) pharm_pert_sens_analysis/run_fig2.m: visualize all trajectories (figure 2) [uploaded]
	Calls pharm_pert_sens_analysis/plot_simulated_traj_fig2.m

(2) PharmPert_CPD.py, PharmPert_FPCA.py, calculate_signaling_codons.m: CPD, fPCA, signaling codon representation calculations
	*For calculate_signaling_codons.m, the updated "oscillation power" calculation is included in get_sig_stats_v2023.m, it 						
	should be added to the NFkB_Codon path once it is downloaded from __

(2.5) organize_CPD_feature_vectors.m, organize_fPCA_feature_vectors.m, organize_codons_feature_vectors.m, organize_traj_cluster.m:
calculate the corresponding vectors for each trajectories (saved as 3000 x num_features matrices) and rearrange vectors into landscapes (heatmaps) to visualize 15 stimulus conditions per drug regime

(3) visualize_landscape_traj_per_drug_regime_fig3.m: visualize landscape of each vector, given user-specified drug regimes (Figure 3)

(3.5) visualize_CPD_comps_figS1.m, visualize_fPCA_temp_patterns_figS2.m: visualize fpca, cpd, components and temporal patterns 

(4) epsilon_network_clustering/optimize_epsilon_network_clustering.m: epsilon network clustering for reduced dimension vectors, to identify the optimal epsilon (Figure 4, table)
	Calls epsilon_network_clustering/epsilon_clustering.m to apply clustering to each drug regime
	Calls epsilon_network_clustering/min_wrong_set_match.m to calculate misclustering rate based on ground truth clusters
	Also relies on the slmetric_pw function in pwmetric folder which was downloaded from 
	https://www.mathworks.com/matlabcentral/fileexchange/15935-computing-pairwise-distances-and-metrics

(4.5) organize_codons_feature_vectors.m: non-responder calculation
      epsilon_network_clustering/find_non_responders_per_drug_regime.m: converts non-responder indices into their corresponding
      stimulus condition index (1-15) under each drug regime

(4.5') xxx.m: non-responder definition for other methods (fpca, cpd, etc)
(Non-responder for other methods are defined as the cluster having the most overlap with the ground_truth non-responder) 

(5) Fig4_02132024.m : Fig5_stimulus_confusion_maps_0520.m: visualize stimulus confusion map (Figure 4 [Xiaolu] & figure S4 [Xiaolu] & Figure 5F-E) 

(6) Fig5_score_maps_20240520.m : score maps (Figure 5A-D)

(7) fig5_supp_traj_codon.m : example trajectories (Figure S5)

(8) Fig6A_0520.m : calculation and visualization heatmap of numbers of drug regime causing stimuli confusion (Figure 6A [uploaded] & epsilon_network_clustering/count_stim_conf_across_drug_regime_ut_v_t_fig6B.m)

(9) epsilon_network_clustering/treated_vs_untreated_clusters_figS6.m: calculation and visualization stimulus confusion map between treated and untreated conditions, example trajectories (Figure S6)
	Call epsilon_network_clustering/make_cluster_squares_treated_vs_untreated.m to generate figure




To be deleted:

(1) revision code Figure4.m:
To be confirmed with [Emily]:
% ground truth non responder definition need to be confirmed (how this is defined?)
non_responder_gt{1} = [];
non_responder_gt{2} = [4,7];
non_responder_gt{3} = [4,7,10];
non_responder_gt{4} = [4];
non_responder_gt{5} = [1:5,7,8,10:14];

02/03/2025 - yes these are correct! Defined using a peak amplitude < 0.05 threshold


(2) revision manuscripts:
[Emily] Signaling codon non-responder should be clarified in methods
[Xiaolu] other methods non-responder should be clarified in methods 
