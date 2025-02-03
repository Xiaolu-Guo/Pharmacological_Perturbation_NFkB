%% Load component weights (lambda) and scales calculated in check_CPD_scale.m 

load("CPD_scale_7c_40c_120423.mat") %tensor_weights_7c, scale_all_7c

%% Load trajectories (bc_sens_table) and codons (metrics)

load('Pharm_metrics_2023osc_051423.mat')
%% Load 7 Component CPD results, generate heatmap and cluster of feature vectors

filename = 'tensor_traj_10drug_norm_7c_noun_allt.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
trajectories = double(tensor_factors{1, 5});

ligands = double(tensor_factors{1, 1});
ligand_doses = double(tensor_factors{1, 2});
drugs = double(tensor_factors{1, 3});
drug_doses = double(tensor_factors{1, 4});

cp_heatmap_7c = make_feature_vector_heatmap(ligands, ligand_doses, drugs, drug_doses, scale_all_7c, tensor_weights_7c);
cp_heatmap_labels = create_heatmap_labels(bc_sens_table);
cp_cluster_7c = reshape_heatmap_to_cluster(cp_heatmap_7c);
cp_cluster_labels = create_cluster_labels(bc_sens_table);


save("CPD_7component_cluster_scaled_0.87VE.mat", "cp_cluster_7c", "cp_cluster_labels", "bc_sens_table", "metrics")
save("CPD_7component_heatmap_scaled_0.87VE.mat", "cp_heatmap_7c", "cp_heatmap_labels", "bc_sens_table", "metrics")

clear sig_tensor_factors
clear tensor_factors
clear trajectories
clear ligands
clear ligand_doses
clear drugs
clear drug doses
%% Load 40 component CPD results, generate heatmap and cluster of feature vectors

filename = 'tensor_traj_10drug_norm_0.95VE_noun_allt_new.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
trajectories = double(tensor_factors{1, 5});
tensor_weights = double(sig_fact_tensor.weights);

ligands = double(tensor_factors{1, 1});
ligand_doses = double(tensor_factors{1, 2});
drugs = double(tensor_factors{1, 3});
drug_doses = double(tensor_factors{1, 4});

cp_heatmap_40c = make_feature_vector_heatmap(ligands, ligand_doses, drugs, drug_doses, scale_all_7c, tensor_weights_7c);
cp_heatmap_labels = create_heatmap_labels(bc_sens_table);
cp_cluster_40c = reshape_heatmap_to_cluster(cp_heatmap_40c);
cp_cluster_labels = create_cluster_labels(bc_sens_table);


save("CPD_40component_cluster_scaled_0.95VE.mat", "cp_cluster_40c", "cp_cluster_labels", "bc_sens_table", "metrics")
save("CPD_40component_heatmap_scaled_0.95VE.mat", "cp_heatmap_40c", "cp_heatmap_labels", "bc_sens_table", "metrics")
%% Calculate scaled/weighted (lambda) CPD feature vectors and organize them into a 150 x (# of drug doses * # of components/features) heatmap

function cp_heatmap = make_feature_vector_heatmap(ligands, ligand_doses, drugs, drug_doses, comp_scales, comp_weights)

    [doses, comps] = size(drug_doses);

    d = 1;
    for i_drugs = 1:10
        for i_ligands = 1:5
            for i_ligand_doses = 1:3
            
                drug_ligand_ligandd = ligands(i_ligands, :) .* ligand_doses(i_ligand_doses, :) .* drugs(i_drugs, :);
            
                drug_ligand_ligandd_scaled = drug_ligand_ligandd .* comp_scales .* comp_weights;
            
                mult_drugd = reshape((drug_ligand_ligandd_scaled .* drug_doses)', [1, doses * comps]); 
            
                cp_heatmap(d, :) = mult_drugd;
            
                d = d + 1;
            
                if d > size(ligands, 1) * size(ligand_doses, 1) * size(drug_doses, 1)
                    break
                end
            end
        end
    end
end

%% Reshape feature vector heatmaps to be compatible with epsilon network clustering (3000 x # of components/features)

function cp_cluster = reshape_heatmap_to_cluster(cp_heatmap)

    i_start = 1;
    i_end = 15;

    c_start = 1;
    c_end = 300;

    ii = 1;
    while ii <= 10
        cp_hm_1drug = cp_heatmap(i_start:i_end,:);
    
        j_start = 1;
        j_end = 7;

        jj_start = 1;
        jj_end = 15;
    
        jj = 1;
        while jj <= 20
            cp_hm_reshape(jj_start:jj_end, :) = cp_hm_1drug(:, j_start:j_end);
        
            j_start = j_start + 7;
            j_end = j_end + 7;

            jj_start = jj_start + 15;
            jj_end = jj_end + 15;
        
            jj = jj + 1;
        end
        cp_cluster(c_start:c_end, :) = cp_hm_reshape;
    
        i_start = i_start + 15;
        i_end = i_end + 15;

        c_start = c_start + 300;
        c_end = c_end + 300;
    
        ii = ii + 1;
        clear cp_hm_1drug
        clear cp_hm_reshape
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