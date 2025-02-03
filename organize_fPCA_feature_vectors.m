%% Load trajectories (bc_sens_table) and codons (metrics)

load('Pharm_metrics_2023osc_051423.mat')
%% Load 5 component fPCA results, generate heatmap and cluster of feature vectors

filename = 'fpca_5c_090623.pkl';
fid = py.open(filename,'rb');
fpca_load = py.pickle.load(fid);
fpca_scores = double(fpca_load.data);

fpca_heatmap = feature_vector_heatmap(fpca_scores);

%fpca_scores is already representative of each condition's feature vector
fpca_cluster = fpca_scores;
fpca_heatmap_labels = create_heatmap_labels(bc_sens_table);
fpca_cluster_labels = create_cluster_labels(bc_sens_table);

save("fPCA_5component_cluster.mat", "fpca_cluster", "fpca_cluster_labels", "bc_sens_table", "metrics")
save("fPCA_5component_heatmap.mat", "fpca_heatmap", "fpca_heatmap_labels", "bc_sens_table", "metrics")
%% Organize feature vectors (i.e., component scores) into a 150 x (# of components * 20 drug doses) heatmap

function fpca_heatmap = feature_vector_heatmap(fpca_scores)
    d = 1;
    [num_conditions, num_comps] = size(fpca_scores);
    og_start = 1;
    og_end = 15;
    hm_row_start = 1;
    hm_row_end = 15;

    for i_drugs = 1:10
        hm_col_start = 1;
        hm_col_end = num_comps;
        for i_ddoses = 1:20
            

            
            fpca_heatmap(hm_row_start:hm_row_end, hm_col_start:hm_col_end)= ...
                        fpca_scores(og_start:og_end, :);
            
            og_start = og_start + 15;
            og_end = og_end + 15;
            
            hm_col_start = hm_col_start + num_comps;
            hm_col_end = hm_col_end + num_comps;
        end
        hm_row_start = hm_row_start + 15;
        hm_row_end = hm_row_end + 15;
    end

end

%% Generate labels for each condition (n = 3000)
function labels = create_cluster_labels(bc_sens_table)
    i_species_names = 2; %NFkBn
    dose_names = unique(bc_sens_table.dose_str);
    bc_names = unique(bc_sens_table.bc_group);
    bc_names_mod = bc_names([10, 1:9]);
    dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
    ligand_names = unique(bc_sens_table.ligand, 'stable');
    species_names = unique(bc_sens_table.species);

    DD = (1:20);
    for i_labels = DD
        DDx{i_labels} = strcat('DD', num2str(DD(i_labels)));
    end

    ligand_names_space = {'-TNF-', '-LPS-', '-CpG-', '-PolyIC-', '-Pam3CSK-'};
    dose_names_space = {'LSD-', 'MSD-', 'HSD-'};

    i_cg = 1;
    while i_cg <= 3000
    for i_drug = 1:length(bc_names_mod)
        for i_dd = 1:length(DDx)
            for i_ligand = 1:length(ligand_names_space)
                for i_ld = 1:length(dose_names_space)
                    cg_label{i_cg} = strcat(bc_names_mod{i_drug}, ligand_names_space{i_ligand},...
                    dose_names_space{i_ld}, DDx{i_dd});
                    i_cg = i_cg + 1;
                end
            end
        end
    end
    end

    labels = cg_label';
end

%% Generate labels for each condition in heatmap format in order to extract 
% all 15 stimulus conditions within a given drug regime

function heatmap_labels = create_heatmap_labels(bc_sens_table)
    i_species_names = 2; %NFkBn
    dose_names = unique(bc_sens_table.dose_str);
    bc_names = unique(bc_sens_table.bc_group);
    bc_names_mod = bc_names([10, 1:9]);
    dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
    ligand_names = unique(bc_sens_table.ligand, 'stable');
    species_names = unique(bc_sens_table.species);

    ligand_names_space = {'-TNF-', '-LPS-', '-CpG-', '-PolyIC-', '-Pam3CSK-'};
    dose_names_space = {'LSD-', 'MSD-', 'HSD-'};
    
    %cluster_labels_hm = zeros(150, 20);
    i_cg = 1;
    while i_cg <= 150
        for i_drug = 1:length(bc_names_mod)
            for i_ligand = 1:length(ligand_names_space)
                for i_ld = 1:length(dose_names_space)
                    hm_label_nodd{i_cg,1} = strcat(bc_names_mod{i_drug}, ligand_names_space{i_ligand},...
                        dose_names_space{i_ld});
                     i_cg = i_cg + 1;
                end
            end
        end
    end

    DD = (1:20);
    for i_repeat = 1:150
        for i_labels = DD
        DDx_hm_label{i_repeat, i_labels} = strcat('DD', num2str(DD(i_labels)));
        end
    end

    for i_col = 1:20
        heatmap_labels(:, i_col) = strcat(hm_label_nodd, DDx_hm_label(:, i_col));
    end

end