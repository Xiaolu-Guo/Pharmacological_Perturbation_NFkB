
%%

load('optimal_clusters_all_methods_CPD_scaled_022724.mat')

load('codon_clusters_0520_reorder.mat')


methods_struct.codons.cluster_array_all
cluster_array_all = methods_struct.codons.cluster_array_all([21:length(cluster_array_all),1:20]);

%%

mycolor_map = [1,1,1;
    1,0,0;
    0,0,0];

% Create a figure

pic_num_each_column =10;
pic_num_each_row = 20;


paperpos=[0,0,100*pic_num_each_row,80*pic_num_each_column]*2;
papersize=[100*pic_num_each_row 80*pic_num_each_column]*2;
ptcolor=[0,0,0];

% figure;

%% change the sub figure size here
subfig_left = 0.1;
subfig_bottom = 0.1;
subfig_width = 0.8;
subfig_length = 0.8;

% figure(1)
% set(gcf, 'PaperUnits','points')
% set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
%
% figure(2)
% set(gcf, 'PaperUnits','points')
% set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)

i_map = 1;
inhibition_score = [];
specificity_score = [];
confusion_score = [];

% reorder the cluster_array
cluster_array_all = methods_struct.codons.cluster_array_all([21:length(cluster_array_all),1:20]);

for i_row = 1:pic_num_each_column% ;%len6

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_column = 1:pic_num_each_row

        b = cluster_array_all{i_map};
        non_responder_cluster_elements = non_responders_array{i_map}';
        inhibition_score(i_row,i_column) = length(non_responder_cluster_elements);

        indexToRemove = find(cellfun(@(x) isequal(x, non_responder_cluster_elements), b));

        % Remove the element at the found index
        if ~isempty(indexToRemove) % Check if the element was found
            b(indexToRemove) = []; % Remove the element
        end

        specificity_score(i_row,i_column) = length(b);
        if isempty(b)
            confusion_score(i_row,i_column) = 0;
        else
            confusion_score(i_row,i_column) = max(cellfun(@length, b));
        end

        i_map = i_map +1;


    end
end

if 1
    paperpos = [0,0,400,200];
    papersize = [400,200];
    figure(4)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
    specificity_score_plot = specificity_score(:,:);%[2:10,1]
    h = heatmap(specificity_score_plot)
    h.ColorbarVisible = 'off'; % Remove colorbar
    caxis([0,15])
    h.XDisplayLabels = repmat({''}, 1, size(specificity_score_plot,2)); % Remove x-axis labels
    h.YDisplayLabels = repmat({''}, 1, size(specificity_score_plot,1));
    saveas(gcf,'../subfigures/specificity_score_0520','epsc')
    close()

    figure(4)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
    confusion_score_plot = confusion_score(:,:);%[2:10,1]
    h = heatmap(confusion_score_plot)
    h.ColorbarVisible = 'off'; % Remove colorbar
    caxis([0,15])
    h.XDisplayLabels = repmat({''}, 1, size(confusion_score_plot,2)); % Remove x-axis labels
    h.YDisplayLabels = repmat({''}, 1, size(confusion_score_plot,1));
    saveas(gcf,'../subfigures/confusion_score_0520','epsc')
    close()

    figure(4)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
    inhibition_score_plot = inhibition_score(:,:);%%[2:10,1]
    h = heatmap(inhibition_score_plot)
    h.ColorbarVisible = 'off'; % Remove colorbar
    caxis([0,15])
    h.XDisplayLabels = repmat({''}, 1, size(inhibition_score_plot,2)); % Remove x-axis labels
    h.YDisplayLabels = repmat({''}, 1, size(inhibition_score_plot,1));
    saveas(gcf,strcat(savefig_path,'inhibition_score_0520'),'epsc')
    close()

end

if 1
    cluster_category(1) = sum(inhibition_score == 15,'all')/length(inhibition_score(:));
    cluster_category(2) = sum(inhibition_score <= 1,'all')/length(inhibition_score(:));
    cluster_category(3) = 1 - cluster_category(1) - cluster_category(2);

    paperpos = [0,0,100,100];
    papersize = [100,100];

    figure(4)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
    ylabel({'percentage of','drug regimes'})
    bar(cluster_category*100)
    ytickformat('percentage')
    set(gca,'XTick',{});
    saveas(gcf,strcat(savefig_path,'cluster_category_counts'),'epsc')
    close()

end
