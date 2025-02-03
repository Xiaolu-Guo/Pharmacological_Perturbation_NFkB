%% Load trajectories (bc_sens_table) and codons (metrics)

load('Pharm_metrics_2023osc_051423.mat')
% metrics = 10 (drugs) x 5 (stimuli) cell
    % inside of each is a 1 x 3 cell representing the 3 stim doses
        % inside of each is a struct consisting of 21 x n vectors
        % calculated for each metric

metric_fields = {'duration', 'oscpower', 'max_value', 'max_pos_integral',  'time2HalfMaxPosIntegral', 'pk1_time'};%EvsL - smaller time, earlier activity
metric_titles = {'Duration', 'Osc power', 'Max value', 'Total activity', 'Early vs late', 'Speed'};

bc_names = unique(bc_sens_table.bc_group);
bc_names_mod = bc_names([10, 1:9]);
ligand_names = unique(bc_sens_table.ligand, 'stable');
mparams = 10.^linspace(0, 3, 21);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};

[metric_cluster, untreated_metric_cluster] = extract_codon_cluster(metrics, bc_names, ligand_names, metric_fields);

% Non-responders are defined as conditions with trajectories that have a
% peak conc. < 0.05 nM
[no_response_row, no_response_col] = find(metric_cluster(:, 3) < 0.05);
metric_cluster(no_response_row, 5) = 4; %Ensuring even EvL activity
metric_cluster(no_response_row, 6) = 0; %Speed = 0

%Z-score codons so they are on the same scale
metric_cluster_norm(:, 1:5) = zscore(metric_cluster(:, 1:5));
total_mean(1:5) = mean(metric_cluster(:, 1:5));
total_std(1:5) = std(metric_cluster(:, 1:5));
total_mean(6) = nanmean(metric_cluster(:, 6));
total_std(6) = nanstd(metric_cluster(:, 6));
metric_cluster_norm(:, 6) = (metric_cluster(:,6) - nanmean(metric_cluster(:, 6))) ...
                           /nanstd(metric_cluster(:, 6));
metric_cluster_norm(isnan(metric_cluster_norm)) = 0;

for i = 1:size(untreated_metric_cluster, 1)
    if untreated_metric_cluster(i, 3) < 0.05
        untreated_metric_cluster(i, 5) = 4; %EvL
        untreated_metric_cluster(i, 6) = 0; %1st peak time
    end
end

untreated_metric_cluster_norm = (untreated_metric_cluster - total_mean) ./ total_std;


cluster_labels = create_cluster_labels(bc_sens_table);

metric_heatmap = reshape_cluster_to_heatmap(metric_cluster_norm);
heatmap_labels = create_heatmap_labels(bc_sens_table);

save("Codon_cluster_norm_nonres_correction.mat", "metric_cluster_norm", "untreated_metric_cluster_norm", ...
    "cluster_labels", "bc_sens_table", "metrics")
save("Codon_heatmap_norm_nonres_correction.mat", "metric_heatmap", "heatmap_labels", "bc_sens_table", "metrics")

%Saving non-responder indices
save("non_responders_inds.mat", "no_response_row", "cluster_labels")
%% Extracts specified n metrics (codons) and organizes them into a 3000 x n cluster

function [metric_cluster, untreated_metric_cluster] = extract_codon_cluster(metrics, bc_names, ligand_names, metric_fields)
c = 1;
while c <= 3000
    for i_bc_group = [10, 1:9] %TSA first
        
        if strcmp(bc_names(i_bc_group), 'TSA')
            i_mparam_range = 2:21;
        else
            i_mparam_range = flip(1:20);
        end
        
        for i_mparam = i_mparam_range

            ut = 1;
            for i_ligand = 1:length(ligand_names)
            
                for i_dose = 1:3
                        for i_metric = 1:length(metric_fields)
                            if metric_fields{i_metric} == "duration"
                                col = 2;
                            else
                                col = 1;
                            end
                            if metric_fields{i_metric} == "time2HalfMaxPosIntegral"
                                m = abs(metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(i_mparam,col) - 8);
                                untreated_m = abs(metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(untreated_ind,col) - 8);
                                
                                %if metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{3})(i_mparam,col) < 0.05
                                %    m = 4;
                                %    untreated_m = 4;
                                %else
                                %    m = abs(metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(i_mparam,col) - 8);
                                %    untreated_m = abs(metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(untreated_ind,col) - 8);
                                %end
                            else
                                m = metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(i_mparam,col);
                                untreated_m = metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(untreated_ind,col);
                            end
                            
                            if metric_fields{i_metric} == "pk1_time" %Speed of first detectable peak
                                if metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{3})(i_mparam,col) < 0.05 
                                    m = 0;
                     
                                else
                                    m = metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(i_mparam,col);
                                    untreated_m = metrics{i_bc_group, i_ligand}{i_dose}.(metric_fields{i_metric})(untreated_ind,col);
                                end
                            end
                                    

        
                            metric_cluster(c, i_metric) = m;
                            untreated_metric_cluster(ut, i_metric) = untreated_m;
                        

                        end
                        ut = ut + 1;
                        c = c + 1;
                 end

            end
        end
    end
end
end

%% reformatting cluster matrix into a 150 x (# of features * 20 drug doses) heatmap

function metric_heatmap = reshape_cluster_to_heatmap(metric_cluster_norm)

    d = 1;
    [num_conditions, num_comps] = size(metric_cluster_norm);
    og_start = 1;
    og_end = 15;
    hm_row_start = 1;
    hm_row_end = 15;

    for i_drugs = 1:10
        hm_col_start = 1;
        hm_col_end = num_comps;
        for i_ddoses = 1:20
            

            
            metric_heatmap(hm_row_start:hm_row_end, hm_col_start:hm_col_end)= ...
                        metric_cluster_norm(og_start:og_end, :);
            
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