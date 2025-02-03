fig_path = '../Fig3/';

%% Load trajectories (bc_sens_table) and codons (metrics)

load('Pharm_metrics_2023osc_051423.mat')

%% Load all feature vector heatmaps
load("CPD_7component_heatmap_scaled_0.87VE.mat") %cp_heatmap_40c
load("CPD_40component_heatmap_scaled_0.95VE.mat") %cp_heatmap_40c
load("fPCA_5component_heatmap.mat") %fpca_heatmap
load("Codon_heatmap_norm_nonres_correction.mat") %metric_heatmap

%%
drug_names = {'TSA', 'BTRCPi', 'CHL', 'CHX', 'IKKi', 'MG132', 'PDTC',...
              'PP2Ai', 'Sel', 'TAK1i'};


%For example drug regimes - Sel @ DD15 and MG132 @ DD6
feature_space = "Codons";
drug_choice_all = {'Sel', 'MG132'};
drug_dose_choice_all = [15, 6]; 
num_features = 6; %7C CPD = 7; 40C CPD = 40; 5C fPCA = 5; codons = 6

drug_regime_heatmap(feature_space, fig_path, metric_cluster, num_features, drug_label,... 
                    drug_choice_all, drug_dose_choice_all);
                
plot_traj_per_drug_regime(bc_sens_table, drug_choice_all(1), drug_dose_choice_all(1), fig_path);
plot_traj_per_drug_regime(bc_sens_table, drug_choice_all(2), drug_dose_choice_all(2), fig_path);
%% Plot all 15 stim condition trajectories (option to include untreated) for a given drug regime
function plot_traj_per_drug_regime(bc_sens_table, bc_choice, dd_choice, fig_path)
figure()
i_species_names = 2; %NFkBn
dose_names = unique(bc_sens_table.dose_str);
bc_names = unique(bc_sens_table.bc_group);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
ligand_names = unique(bc_sens_table.ligand, 'stable');
species_names = unique(bc_sens_table.species);


%ligand_color = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880],...
%                [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980]};

ligand_color = {[119 180 202]/255, [222 78 66]/255, [137 180 66]/255, ...
    [101 77 123]/255, [229 129 56]/255};

pic_num_each_row = 3;
pic_num_each_column = 5;
subfig_left = 0.10;
subfig_bottom = 0.05;
subfig_width = 0.80;
subfig_length = 0.70;

for i_bc_names = 1:length(bc_choice)
    figure()

        for i_ligand_names = 1:length(ligand_names)
            
            d = 1;
            for i_dose_names = 1:length(dose_names)
               
                index = find(strcmp(bc_sens_table.bc_group,bc_choice{i_bc_names}) &...
                    strcmp(bc_sens_table.dose_str, dose_names{i_dose_names}) &...
                    strcmp(bc_sens_table.ligand, ligand_names{i_ligand_names}) &...
                    strcmp(bc_sens_table.species,species_names{i_species_names}));
                if bc_choice == "TSA"
                    untreated_ind = 1;
                else
                    untreated_ind = 21;
                end
               
                if isempty(index) ~= 1
                    for i_index = 1:length(index)
                        i_sim_data = index(i_index);
                        
                        for i_dd = 1:length(dd_choice)
                        
                        if bc_sens_table.multiplier_count{i_sim_data} == dd_choice(i_dd)
                            
                         position_subfig=[1/pic_num_each_row*(d - 1)+0.5/pic_num_each_row*subfig_left,...
                          (1 - 0.9/pic_num_each_column * i_ligand_names),...
                            1.175/pic_num_each_row*subfig_width,...
                            1.15/pic_num_each_column*subfig_length];
                        
                            subplot('Position',position_subfig)
                            plot((1:length(bc_sens_table.trajectory{i_sim_data}))/60,bc_sens_table.trajectory{i_sim_data},...
                                'LineWidth', 2, 'color', ligand_color{i_ligand_names});
                            hold on 
                            
                            %FOR WHEN WANT UNTREATED INCLUDED IN PLOT
                            plot((1:length(bc_sens_table.trajectory{i_sim_data}))/60,bc_sens_table.trajectory{index(untreated_ind)},'--',...
                                'LineWidth', 1, 'color', 'k');
                            %plot((1:length(bc_sens_table.trajectory{i_sim_data}))/60, trajectories(:, 2));
                            ylim([0, 0.3])
                            xlim([0 8])
                            set(gca, 'xtick', [0, 4, 8], 'xticklabels', [], ...
                                'ytick', [0, 0.15, 0.3], 'yticklabels', [])
                            ax=gca;
                            ax.LineWidth = 2;
                            %ylabel(ligand_names(i_ligand_names), 'Fontweight', 'bold')
                            %xlabel('Time (hrs)')
                            grid on
                           
                            
                        end
                        end
                        
                    end
                    
                     d = d + 1;
                end
            end
        end
end

drug_dd_label = strcat(bc_choice, "DD", num2str(dd_choice)); 
saveas(gcf,strcat(fig_path,"15_stim_", drug_dd_label_, "treated_untreated_traj"),'png')
end
                
%%
function drug_regime_heatmap(feature_space, fig_path, heatmap, num_features, drug_label, drug_choice, drug_dose_choice)
    cc = 1;
    i_start = 1;
    i_end = num_features;
    while cc <= 20
        for ii = i_start:i_end
            dd_label(1, ii) = cc;
        end
        cc = cc + 1;
        i_start = i_start + num_features;
        i_end = i_end + num_features;
    end
    
    %All drug regimes (i.e., drug_choice/drug_dose_choice combos) are
    %plotted on the same heatmap so their scales/color ranges are the same!
    col_start = 1;
    col_end = num_features;
    for ii = 1:length(drug_choice)
        dd_ind = find(dd_label == drug_dose_choice(ii));
        drug_ind = find(strcmp(drug_label, drug_choice(ii)));

        dr_hm(:,col_start:col_end) = heatmap(drug_ind, dd_ind);
        col_start = col_start + num_features;
        col_end = col_end + num_features;
        drug_dd_label(ii) = strcat(drug_choice(ii), "DD", num2str(drug_dose_choice(ii)), "_"); 
    end
    
    figure()
    imagesc(dr_hm)
    hold on
    colormap(hot)
    c = colorbar('TickLabels', [], 'FontWeight', 'Bold'); %'Ticks', [], 'TickLabels', [],
    set(gca, 'TickLength', [0 0],'ytick', [], ...
        'yticklabel', [], 'xtick', [], 'xticklabel', [])
    
    drug_dd_label_combined = {cat(2, drug_dd_label{:})};
    
    saveas(gcf,strcat(fig_path,feature_space, "_space_", drug_dd_label_combined, "landscapes"),'png')
end                