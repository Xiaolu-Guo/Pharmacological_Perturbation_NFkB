%% Load 5 component fPCA results
filename = '../fpca_5c_components_090623.pkl';
fid = py.open(filename,'rb');
fpca_load = py.pickle.load(fid);
fpca_components = double(fpca_load.data_matrix);
[num_comps, timepoints] = size(fpca_components);


%% Temporal pattern plots
x = 1:1:timepoints;

traj_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...	
              [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],... 
              [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
figure()
for i = 1:num_comps
    subplot(4, 4, i + 1)
    plot(x, fpca_components(i,:), 'LineWidth', 2, 'color', traj_colors{i});
    %set(gca, 'xtick', [0:240:480], 'xticklabels', [], ...
    %            'ytick', [-4:2:4], 'yticklabels', [])
    ylim([-4, 4])
    xlim([0, 481])
    set(gca, 'xtick', [], 'ytick', [])
    ax=gca;
    ax.LineWidth = 2;
end

fig_path = '../FigS2/';
saveas(gcf,strcat(fig_path,'fPCA_5comp_temp_patterns'),'png')