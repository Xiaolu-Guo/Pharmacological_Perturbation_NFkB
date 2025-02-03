%% 
% Path for all codon calculation functions
codon_path = '../NFkB_codon/';

% Load simulation results
load('Pharm_all_stim_0.001-1X_simulations.mat')

dose_names = unique(bc_sens_table.dose_str);
bc_names = unique(bc_sens_table.bc_group);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
ligand_names = unique(bc_sens_table.ligand, 'stable');
species_names = unique(bc_sens_table.species);

c = 0;
z = 0;
v = 0;


vis_data_field = {'model_sim'};%,'sample'};
data_label = {'simulation'};%,'sample'};

for i_bc_names = 1:length(bc_names)

     for i_ligand = 1:length(ligand_names)
         ligand_ind = contains(bc_sens_table.ligand, ligand_names{i_ligand});
         bc_sens_table_1l = bc_sens_table(ligand_ind,:);
 
        for i_dose_names = 1:length(dose_names)
         
           for i_species_names = 1:length(species_names)
                 index = find(strcmp(bc_sens_table_1l.bc_group,bc_names{i_bc_names}) &...
                     strcmp(bc_sens_table_1l.dose_str,dose_names{i_dose_names}) &...
                     strcmp(bc_sens_table_1l.species,species_names{i_species_names}));
             
                  i_cell = 1;
             
                 for i_index = 1:length(index)
                     i_sim_data = index(i_index);
                     i_para = bc_sens_table_1l.multiplier_count{i_sim_data};
                 
                     a=bc_sens_table_1l.trajectory{i_sim_data};
                     sim{i_species_names}(i_cell,:)  = a(1:5:end);
%                 %
%                 %                 plot((1:length()/60,bc_sens_table.trajectory{i_sim_data},...
%                 %                     'LineWidth',Line_wid(i_para),'Color',color_mapping(i_para,:));hold on
                     i_cell = i_cell + 1;
%                 
                 end
             
             end
 
             if ~isempty(index)
                c = c + 1;
%                 % NFkBn + IkBaNFkBn
                 data.model_sim{c} =sim{2};
                 data.info_ligand{c} = bc_sens_table_1l.ligand{i_sim_data};
                 data.info_dose_index{c} = i_dose_names;
                 data.info_dose_str{c} = bc_sens_table_1l.dose_str{i_sim_data};
                 data.info_num_cells{c} = bc_sens_table_1l.multiplier_count{i_sim_data};
 
                 data.order{c} = (1:data.info_num_cells{c})';
                 data.time_points_num = 97;
             end
        end
         c = 0;
% 
%     
% 
         data.exp = data.model_sim;
%     
         cd(codon_path)
% 
         [collect_feature_vects{i_bc_names, i_ligand},metrics{i_bc_names, i_ligand}] = calculate_codon(data,vis_data_field,data_label);%,  parameter
%         
         cd(currentfolder)
% 
     end
end

save("Pharm_metrics_2023osc_051423.mat", "metrics", "bc_sens_table")