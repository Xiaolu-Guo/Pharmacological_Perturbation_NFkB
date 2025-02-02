fig_path = '../FigS1/';

%% Load 7 component CPD results
filename = '../tensor_traj_10drug_norm_7c_noun_allt.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
trajectories = double(tensor_factors{1, 5});
stim_comps = double(tensor_factors{1, 1});
dstim_comps = double(tensor_factors{1, 2});
drug_comps = double(tensor_factors{1, 3});
ddrug_comps = double(tensor_factors{1, 4});
%% component figures

%Colors chosen to highlight high dose TNF @ TSA drug dose 2 condition's
%weightings in component 6. Can be updated depending on condition/component
%of interest

for cl = 1:7
    comp_legend{cl} = strcat('Comp.', num2str(cl));
end

dimension_names = ["Stimulus", "Stimulus_dose", "Drug", "Drug_dose", "Temporal_patterns"];

for jj = 1:5
    figure()
    comps = double(tensor_factors{1, jj});

    width = 0.1;
    width_plus = -3:1:3;
    %subplot(1, 4, jj)
    colormap(hot)
    x = 1:1:size(comps, 1);
    for ii = 1:size(stim_comps, 2)
        if jj == 4
            if ii == 6
                plot(x, comps(:, ii), 'LineWidth', 2.5, 'color', [0.9290 0.6940 0.1250])
                %plot(5, comps(5, 3), 'k*', 'Linewidth', 3)
                hold on
            else
                plot(x, comps(:, ii), 'LineWidth', 2.5, 'color', 'k')
                hold on
            end
            set(gca, 'xtick', [1:1:20], 'xticklabels', [], ...
                 'ytick', [0:0.2:1.0], 'yticklabels', [])
            [lgnd,icons,plots,txt] = legend(comp_legend);
            temp = [lgnd; lgnd.ItemText];
            set(temp,'FontSize',10);
            set(temp,'FontName','Calibri', 'FontWeight', 'Bold');
            ylim([0, 1.0])
        elseif jj == 5
            plot(x, comps(:, ii), 'LineWidth', 2.5)
            hold on
            set(gca, 'xtick', [0:240:480], 'xticklabels', [], ...
                'ytick', [0:0.1:0.3], 'yticklabels', [])
            ylim([0, 0.3])
        else
            if ii == 3
                bar((x + (width_plus(ii) * width)), comps(:, ii), width, 'FaceColor',[0.9290 0.6940 0.1250])
                hold on
            else
                bar((x + (width_plus(ii) * width)), comps(:, ii), width, 'FaceColor', 'k', 'EdgeColor', [0 0.4470 0.7410], 'Linewidth', 1.5)
                hold on
            end

            set(gca, 'xtick', [], 'ytick', [0:0.2:1.0], 'yticklabels', [])
            ylim([0, 1.0])
        end
        
    end

    xlim([0.5, (size(comps, 1) + 0.5)])
    ax=gca;
    ax.LineWidth = 1.5;
    
    saveas(gcf,strcat(fig_path, dimension_names(jj), '_component_weights_7comp_CPD'),'png')
end

%% Variance explained/error plot
figure()
x = 1:1:8;
var_explained = [0.6358028558021684, 0.748831, 0.7943326, 0.809362207, 0.842756195, 0.85442729, 0.87515641, 0.8869489150111896];
plot(x, var_explained, 'LineWidth', 2, 'color', 'k')
ylim([0, 1.0])
set(gca, 'xtick', [0:1:8], 'xticklabels', [], ...
    'ytick', [0:0.2:1.0], 'yticklabels', [])
ax=gca;
ax.LineWidth = 1.5;
saveas(gcf,strcat(fig_path,'CPD_error_plot_comp1-8'),'png')

%% Randomized timepoint order temporal pattern plots

filename = '../tensor_random_timepoint_order.pkl';
fid = py.open(filename,'rb');
sig_fact_tensor = py.pickle.load(fid);
tensor_factors = cell(sig_fact_tensor.factors);
tensor_weights = double(sig_fact_tensor.weights);
trajectories = double(tensor_factors{1, 5});
stim_comps = double(tensor_factors{1, 1});
dstim_comps = double(tensor_factors{1, 2});
drug_comps = double(tensor_factors{1, 3});
ddrug_comps = double(tensor_factors{1, 4});

x = 1:1:481;
traj_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...	
              [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],... 
              [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};
          
traj_order = [1, 2, 5, 4, 7, 6, 3];          

figure()
for i_traj = 1:size(trajectories, 2)
    
    subplot(4, 4, i_traj + 1)
  
    plot(x, trajectories(:, traj_order(i_traj)), 'LineWidth', 2, 'color', traj_colors{i_traj})
    ylim([0, 0.2])
    xlim([0, 481])
    set(gca, 'xtick', [], 'ytick', [])
    ax=gca;
    ax.LineWidth = 2;
    
    %hold on
end

set(gcf, 'PaperUnits', 'inches'); 
set(gcf, 'PaperSize', [8 8]);

saveas(gcf,strcat(fig_path,'CPD_rand_timepoint_order_temp_patterns'),'png')