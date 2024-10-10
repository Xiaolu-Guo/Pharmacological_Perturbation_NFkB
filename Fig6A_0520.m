
%
% a = zeros(15,15);
% a(2,1) = 1;
% a(7,5) = 1;
% a(12,4) = 0.5;
% mycolor_map = [1,1,1;
%                0.7,0.7,0.7;
%                0,0,0];
%
% h = heatmap(a, 'ColorMap', mycolor_map);
% h.ColorbarVisible = 'off'; % Remove colorbar
% h.XDisplayLabels = repmat({''}, 1, size(a,2)); % Remove x-axis labels
% h.YDisplayLabels = repmat({''}, 1, size(a,1)); % Remove y-axis labels
% h.CellLabelColor = 'none'; % Remove text/numbers inside cells
% caxis([0,1])
% grid off

%%
%load('codon_clusters_011924.mat')
% load('optimal_clusters_all_methods_CPD_scaled_120623.mat')
% cluster_array_all = cluster_array_all([21:length(cluster_array_all),1:20]);
% non_responders_array = non_responders_array([21:length(non_responders_array),1:20]);
% treatment_label = treatment_label([(20*15+1):length(treatment_label), 1:(20*15)]);
% save('codon_clusters_0520_reorder.mat','cluster_array_all','non_responders_array','treatment_label');

load('codon_clusters_0520_reorder.mat')

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


%% change the sub figure size here
subfig_left = 0.1;
subfig_bottom = 0.1;
subfig_width = 0.8;
subfig_length = 0.8;

i_map = 1;
A_sum = zeros(15);
for file_num=1:pic_num_each_column% ;%len6
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:pic_num_each_row
        
        b = cluster_array_all{i_map};
        non_responder_cluster_elements = non_responders_array{i_map};
        
        % Initialize matrix a
        a = zeros(15, 15);
        
        % Mark diagonal elements (each element is always clustered with itself)
        a(1:15+1:end) = 1;
        
        a_full = a;
        
        % Iterate over each cell in b
        for i = 1:length(b)
            % If the cell contains more than one element
            if numel(b{i}) > 1
                % Get all pairs of elements within this cluster
                cluster_elements = b{i};
                pairs = nchoosek(cluster_elements, 2); % All combinations of 2 elements
                
                % Mark each pair as clustered in matrix a
                for j = 1:size(pairs, 1)
                    if pairs(j, 2)>pairs(j, 1)
                        a(pairs(j, 2), pairs(j, 1)) = 1;
                    else
                        a(pairs(j, 1), pairs(j, 2)) = 1;
                    end
                    a_full(pairs(j, 2), pairs(j, 1)) = 1;
                    a_full(pairs(j, 1), pairs(j, 2)) = 1;
                    
                end
            end
        end
        
        % non responder
        for j = 1:length(non_responder_cluster_elements)
            a(non_responder_cluster_elements(j),non_responder_cluster_elements(j)) = 0.5;
            a_full(non_responder_cluster_elements(j),non_responder_cluster_elements(j)) = 0.5;
            
        end
        
        if length(non_responder_cluster_elements)>1
            pairs_non_responser = nchoosek(non_responder_cluster_elements, 2); % All combinations of 2 elements
            
            for j = 1:size(pairs_non_responser, 1)
                if pairs_non_responser(j, 2)>pairs_non_responser(j, 1)
                    a(pairs_non_responser(j, 2), pairs_non_responser(j, 1)) = 0;
                else
                    a(pairs_non_responser(j, 1), pairs_non_responser(j, 2)) = 0;
                    
                end
                a_full(pairs_non_responser(j, 2), pairs_non_responser(j, 1)) = 0;
                a_full(pairs_non_responser(j, 1), pairs_non_responser(j, 2)) = 0;
                
            end
        end
        
        % Create heatmap
        a(triu(true(size(a)), 1)) = NaN;
        A_sum = A_sum + (a == 1);
        
        
        i_map = i_map+1;
    end
end

paperpos = [0,0,100,100];
papersize = [100,100];

figure(4)
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
for i_diag = 1:15
    A_sum(i_diag,i_diag) = 0;
end
imagesc(A_sum);
colormap(parula); % Apply your color map
caxis([0,100]);
colorbar;

%         h = heatmap(a, 'ColorMap', mycolor_map);
%         h.ColorbarVisible = 'off'; % Remove colorbar
%         h.XDisplayLabels = repmat({''}, 1, size(a,2)); % Remove x-axis labels
%         h.YDisplayLabels = repmat({''}, 1, size(a,1)); % Remove y-axis labels
%         h.CellLabelColor = 'none'; % Remove text/numbers inside cells
%         caxis([0,1]);
%         grid off;
%  axis off;

set(gca, 'box', 'off'); % Turn off the box surrounding the plot
%         set(gca, 'XAxisLocation', 'bottom'); % Keep the x-axis (bottom)
%         set(gca, 'YAxisLocation', 'left'); % Keep the y-axis (left)
set(gca, 'XTick', []); % Optionally remove x-axis ticks if desired
set(gca, 'YTick', []); % Optionally remove y-axis ticks if desired


figure(4)
exportgraphics(gcf, '../subfigures/Fig6A_0520.eps','Resolution',600)
% saveas(gcf,'../subfigures/stimuli_confusion_map','eps');
close

