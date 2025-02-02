fig_path = '../FigS2/';
%% Load 5 component fPCA results
filename = 'fpca_5c_components_090623.pkl';
fid = py.open(filename,'rb');
fpca_load = py.pickle.load(fid);
fpca_components = double(fpca_load.data_matrix);
[num_comps, timepoints] = size(fpca_components);


%% Component plots
x = 1:1:timepoints;
figure()
for i = 1:num_comps
    plot(x, fpca_components(i,:), 'LineWidth', 2.5);
    set(gca, 'xtick', [0:240:480], 'xticklabels', [], ...
                'ytick', [-4:2:4], 'yticklabels', [])
    ylim([-4, 4])
    xlim([0, 481])
    hold on 
end
ax=gca;
ax.LineWidth = 1.5;

for i = 1:num_comps
    comp_legend{i} = strcat('Comp.', num2str(i));
end
%[lgnd,icons,plots,txt] = legend(comp_legend);
%temp = [lgnd; lgnd.ItemText];
%set(temp,'FontSize',10);
%set(temp,'FontName','Arial');

saveas(gcf,strcat(fig_path,'fPCA_5comp_temp_patterns'),'png')
%% Variance explained/error plot
filename = 'fpca_var_explained_092023.pkl';
fid = py.open(filename,'rb');
fpca_load = py.pickle.load(fid);
fpca_ve = double(fpca_load.data);

figure()
x = 1:1:length(fpca_ve(1:10));
plot(x, fpca_ve(1:10), '-o', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'color', 'k')
ylim([0, 1.2])
set(gca, 'xtick', [1:length(fpca_ve(1:10))], 'xticklabels', [], ...
    'ytick', [0:0.2:1.2], 'yticklabels', [])
ax=gca;
ax.LineWidth = 1.5;

saveas(gcf,strcat(fig_path,'fPCA_error_plot'),'png')