function [] = plot_simulated_traj_fig1(data_info,figure_info)
addpath('../CommonUsedFunction/')


% draw the highest dose

if nargin <1
    Module_sens = 'TNF';
end


if exist(figure_info.save_figure_path,'dir')
else
    mkdir(figure_info.save_figure_path)
end

load(strcat(data_info.save_file_path,data_info.save_file_name))

%Change color mapping based on biochem vs pharm!
length_para_vec = 21;

%biochem
%color_mapping = zeros(length_para_vec,3);
%color_mapping(1:10,:) = [linspace(1,0.5,floor(length_para_vec/2))',zeros(floor(length_para_vec/2),2)];
%color_mapping(12:21,:) = [zeros(floor(length_para_vec/2),2),linspace(0.5,1,floor(length_para_vec/2))'];

%pharm
%color_mapping_og = zeros(length_para_vec,3);
%color_mapping_og(1:10,:) = [linspace(1,0.5,floor(length_para_vec/2))',zeros(floor(length_para_vec/2),2)];
%color_mapping_og(11:20,:) = [linspace(0.5,0,floor(length_para_vec/2))', zeros(floor(length_para_vec/2),2)];

%color_mapping_tsa = flip(color_mapping_og);

color_mapping_og = flip([ linspace(0, 0.95, 10)',zeros(10,1),zeros(10,1); ...
                    ones(11,1),linspace(0, 0.5, 11)',linspace(0, 0.5, 11)']);
color_mapping_tsa = flip(color_mapping_og);

Line_wid =ones(length_para_vec,1);
Line_wid(floor(length_para_vec/2)+1)=2;

pic_num_each_column = 15;
pic_num_each_row = 5;
%change the sub figure size here
subfig_left = 0.10;
subfig_bottom = 0.05;
subfig_width = 0.85;
subfig_length = 0.70;

paperpos=[0,0,105*pic_num_each_row,80*pic_num_each_column];
papersize=[105*pic_num_each_row 80*pic_num_each_column];
ptcolor=[0,0,0];

bc_names = unique(bc_sens_table.bc_group);
dose_names = {'0.3', '0.33', '3.3', '10', '33', '100', '330'};
ligand_names = unique(bc_sens_table.ligand,'stable');


set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)




if isfield(data_info,'species')

    species_names = data_info.species;
else
    species_names = unique(bc_sens_table.species);
end
i_species_names = 2; %NFkBn


for i_bc_names = 1:length(bc_names)
    if i_bc_names == 1
        figure(1)
        dd = 0;
        d = 0;
        next_bc_set = 0;
    elseif i_bc_names == 6
        figure(2)
        dd = 0;
        d = 0;
        next_bc_set = 0;
    end
   
    
    if i_bc_names == 10
        color_mapping = color_mapping_tsa;
    else
        color_mapping = color_mapping_og;
    end

        for i_ligand_names = 1:length(ligand_names)
            %dd = [0 2 3 4 5];
            next_dose = 0;
            d = d + 1;
            dd = 0;
            for i_dose_names = 1:length(dose_names)
               
                index = find(strcmp(bc_sens_table.bc_group,bc_names{i_bc_names}) &...
                    strcmp(bc_sens_table.dose_str, dose_names{i_dose_names}) &...
                    strcmp(bc_sens_table.ligand, ligand_names{i_ligand_names}) &...
                    strcmp(bc_sens_table.species,species_names{i_species_names}));
               %+ 0.1/pic_num_each_column*subfig_bottom              
                
                if isempty(index) ~= 1
                    next_dose = next_dose + 1;
                    %subplot(15, 5, d + dd);
                    %Set_figure_size
                    %disp(d)
                    position_subfig=[1/pic_num_each_row*(i_ligand_names-1)+0.85/pic_num_each_row*subfig_left,...
                          (1 - 0.9/pic_num_each_column * next_dose - next_bc_set),...
                            1/pic_num_each_row*subfig_width,...
                            1.15/pic_num_each_column*subfig_length];
                        
                    subplot('Position',position_subfig)
                    
                    for i_index = 1:length(index)
                        i_sim_data = index(i_index);
                        i_para = bc_sens_table.multiplier_count{i_sim_data};
                        
                        plot((1:length(bc_sens_table.trajectory{i_sim_data}))/60,bc_sens_table.trajectory{i_sim_data},...
                        'LineWidth',Line_wid(i_para),'Color',color_mapping(i_para,:));hold on
                        
                    end
                
                
                    sti = ligand_names{i_ligand_names};
                    %sti = bc_sens_table.ligand{index(1)};
                    %para_module  = bc_sens_table.parameter_module{index(1)};
                    %ti = title(sti,'FontSize',8);
                    %if any([1:5] == d)
                        %title(ligand_names{i_ligand_names},'FontWeight','Bold','FontSize',16)
                    %end
                    
                    %if any(1:15:61 == d)
                        %ylabel('Low stim. dose','FontWeight','Bold','FontSize',16)
                    %end
                    
                    %if any(6:15:66 == d)
                        %ylabel('Medium stim. dose','FontWeight','Bold','FontSize',16)
                    %end
                    
                    %if any(11:16:71 == d)
                        %ylabel('High stim. dose','FontWeight','Bold','FontSize',16)
                    %end
                        
                    %xlabel('Time(hour)')
                    set(gca,'FontSize',8, 'ytick', [],'xtick', [])
                    xlim([0,8])
                    %text(230,0.27,)
                    if species_names{i_species_names} == "NFkBn"
                        ylim([0,0.3])
                    end
            
                else
                    continue
                end
                dd = dd + 5;
            end

        
        end
        
        next_bc_set = next_bc_set + 1/5;
        
        if i_bc_names == 5
            fig_name = 'Pharm1-5';
            out = char(strcat(figure_info.save_figure_path,fig_name,'_0.001-1X_',species_names{i_species_names}));
            disp(out)
            savefig(gcf,[out '.fig'])

            close
        elseif i_bc_names == 10
            fig_name = 'Pharm6-10';
            out = char(strcat(figure_info.save_figure_path,fig_name,'_0.001-1X_',species_names{i_species_names}));
            disp(out)
            saveas(gcf, out, 'png.')

            close
        end
        d = d + 10;
end
