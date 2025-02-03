sa_choice = 2;
data_ttl = '08-Jun-2022_Pharm_all_stim_0.001-1X_FINAL.mat';
%data_info.species = species2save;
data_info.save_file_path = './data/';
data_info.save_file_name = data_ttl;
if sa_choice == 1
    figure_info.save_figure_path = 'Biochem_figures/NFkB_traj/';
elseif sa_choice == 2
    figure_info.save_figure_path = 'Pharm_figures/NFkB_traj/Fig1/';
end


Fig1(data_info,figure_info);


%%
mparams = 10.^linspace(-3, 0, 21);
for i = 1:length(mparams)
    fcx{i} = strcat(num2str(mparams(i)), 'X');
end

length_para_vec = 21;
color_mapping_og = zeros(length_para_vec,3);
color_mapping_og(1:10,2:3) = [linspace(0.5,0.95,floor(length_para_vec/2))',linspace(0.5,0.95,floor(length_para_vec/2))'];
color_mapping_og(1:10, 1) = ones(10,1);
color_mapping_og(11:20,:) = [linspace(1,0,floor(length_para_vec/2))', zeros(floor(length_para_vec/2),2)];
mymap = flip([ linspace(0, 0.95, 10)',zeros(10,1),zeros(10,1); ones(11,1),linspace(0, 0.5, 11)',linspace(0, 0.5, 11)']);
figure()
colormap(mymap)

set(gca,'ColorScale','log')
caxis([10^-3, 1])
c = colorbar('Ticks', mparams);
DD = (1:20);
for i_labels = 1:length(DD)
    DDx{i_labels+1} = strcat(num2str(DD(i_labels)));
end
DDx{1} = 'Untreated';
c = colorbar('Ticks', mparams, 'TickLabels', []);%flip(DDx),'Fontsize', 8);

