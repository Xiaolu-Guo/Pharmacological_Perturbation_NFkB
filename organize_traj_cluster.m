%% Load trajectories (bc_sens_table) and codons (metrics)

load('Pharm_metrics_2023osc_051423.mat')

%%
bc_names = unique(bc_sens_table.bc_group);
bc_names_mod = bc_names([10, 1:9]);
ligand_names = unique(bc_sens_table.ligand, 'stable');
mparams = 10.^linspace(0, 3, 21);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};

i_species_names = 2; %NFkBn
dose_names = unique(bc_sens_table.dose_str);
%bc_names = unique(bc_sens_table.bc_group);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
ligand_names = unique(bc_sens_table.ligand, 'stable');
species_names = unique(bc_sens_table.species);

% Extract all trajectories, organize them into 3000 x 481 (timepoints)
% matrix for clustering
c = 1;
while c <= 3000
    for i_bc_names = 1:length(bc_names_mod)
    
        if strcmp(bc_names_mod(i_bc_names), 'TSA')
            i_mparam_range = 2:21;
        else
            i_mparam_range = flip(1:20);
        end
        
        DD = 0;
        for i_mparam = i_mparam_range
            DD = DD + 1;

            for i_ligand_names = 1:length(ligand_names)
            
            

                for i_dose_names = 1:length(dose_names)
               
                    index = find(strcmp(bc_sens_table.bc_group,bc_names_mod{i_bc_names}) &...
                        strcmp(bc_sens_table.dose_str, dose_names{i_dose_names}) &...
                        strcmp(bc_sens_table.ligand, ligand_names{i_ligand_names}) &...
                        strcmp(bc_sens_table.species,species_names{i_species_names}));
               
                
                    if isempty(index) ~= 1
                        %if strcmp(bc_names(i_bc_names), 'TSA')
                        %    index = index;
                        %else
                        %    index = flip(index);
                        %end

                        for i_index = 1:length(index)
                            i_sim_data = index(i_index);
                        
                        
                            if bc_sens_table.multiplier_count{i_sim_data} == i_mparam
                        
                                label = convertCharsToStrings(strcat(bc_names_mod{i_bc_names}, '-DD', num2str(DD), '-',ligand_names{i_ligand_names}, '-', dose_names{i_dose_names}));
                                cluster_labels{c, 1} = label;
                                traj_cluster(c, :) = bc_sens_table.trajectory{i_sim_data}';
                                
                                c = c + 1;
                            
                            end
                         end
                        
                     end
                    end
                end

            end
        end
end

c = 1;
i_bc_names = 2; %Random drug to select untreated
untreated_ind = 21;
while c <= 15
    for i_ligand_names = 1:length(ligand_names)
        for i_dose_names = 1:length(dose_names)
               
             index = find(strcmp(bc_sens_table.bc_group, bc_names_mod{i_bc_names}) &...
                          strcmp(bc_sens_table.dose_str, dose_names{i_dose_names}) &...
                          strcmp(bc_sens_table.ligand, ligand_names{i_ligand_names}) &...
                          strcmp(bc_sens_table.species,species_names{i_species_names}));
               
                
                    if isempty(index) ~= 1
                        %if strcmp(bc_names(i_bc_names), 'TSA')
                        %    index = index;
                        %else
                        %    index = flip(index);
                        %end

                        for i_index = 1:length(index)
                            i_sim_data = index(i_index);
                        
                        
                            if bc_sens_table.multiplier_count{i_sim_data} == untreated_ind
                        
                                traj_cluster_untreated(c, :) = bc_sens_table.trajectory{i_sim_data}';
                                
                                c = c + 1;
                            
                            end
                         end
                        
                     end
        end
    end

end

save("Traj_cluster.mat", "traj_cluster", "traj_cluster_untreated", "cluster_labels");