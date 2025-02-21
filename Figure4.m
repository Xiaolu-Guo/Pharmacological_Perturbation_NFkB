%load('optimal_clusters_all_methods_CPD_scaled_02072024_new.mat')
% load('optimal_clusters_all_methods_TA_thresh_021224.mat')
load('optimal_clusters_all_methods_CPD_scaled_20250122.mat')

paperpos=[0,0,100,100]*2;
papersize=[100 100]*2;


mycolor_map = [1,1,1;
    1,0,0;
    0,0,0];

% old version of non responders
% load('non_responders_per_drug_regime_TA_thresh_021224.mat')
%
example_ind = [48, 91, 106, 175, 196];

% new version of non responders: 
load('codon_clusters_0520_reorder.mat')
example_ind_reordered = [28, 71, 86, 155, 176];

ground_truth_array{1} = {[1],[2],[3],[4],[5,6,7,8,9],[10],[11,12], [13,14,15]};%CHL DD8,  no non-responder
ground_truth_array{2} = {[1],[2,5,9,14,15], [3,6], [4,7],[8,13] [10], [11,12]}; % 5/9 2 14/15 % IKKi DD11  [4,7] non-responder
ground_truth_array{3}= {[1,2,3,13,14,15],[4,7,10],[5,8],[6],[9],[11,12]}; % MG132 DD6 [4,7,10] non-responder
ground_truth_array{4} = {[1], [2], [3], [4], [5,8], [6], [7], [9,14,15], [10], [11,12], [13]}; % Sel DD15; [4,7] non-responder; %7 might be non responder, 14 has a big trough... mayb eby itself
ground_truth_array{5} = {[1:5,7,8,10:14], [6], [9,15]}; %TAK1i DD16; non responder [1:5,7,8,10:14]

non_responder_gt{1} = [];
non_responder_gt{2} = [4,7];
non_responder_gt{3} = [4,7,10];
non_responder_gt{4} = [4];
non_responder_gt{5} = [1:5,7,8,10:14];

non_responder.codons{1} = [];
non_responder.comp40{1} = [];
non_responder.comp7{1} = [];
non_responder.traj{1} = [];

%[4,7]
non_responder.codons{2} = [4,7];
non_responder.comp40{2} = [4,7];
non_responder.comp7{2} = [4,7];
non_responder.traj{2} = [4,7];

%[4,7,10]
non_responder.codons{3} = [4,7,10];
non_responder.comp40{3} = [4,7,10];
non_responder.comp7{3} = [1,13,4,7];
non_responder.traj{3} = [4,7];

% [4]
non_responder.codons{4} = [4];
non_responder.comp40{4} = [4,7];
non_responder.comp7{4} = [4,7,10];
non_responder.traj{4} = [4,7];

%[1:5,7,8,10:14]
non_responder.codons{5} = [1:5,7,8,10:14];
non_responder.comp40{5} = [1:15];
non_responder.comp7{5} = [1:15];
non_responder.traj{5} = [1:5,7:14];

%%
if 1
    %%
    
    for i_ex = 1:length(example_ind_reordered)
        clusters_to_comp = ground_truth_array{i_ex};
        
        b = clusters_to_comp;
        non_responder_cluster_elements = non_responder_gt{i_ex};
        
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
        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        
        a(triu(true(size(a)), 1)) = NaN;
        imagesc(a);
        colormap(mycolor_map); % Apply your color map
        caxis([0,1]);
        set(gca, 'box', 'off'); % Turn off the box surrounding the plot
        set(gca, 'XAxisLocation', 'bottom'); % Keep the x-axis (bottom)
        set(gca, 'YAxisLocation', 'left'); % Keep the y-axis (left)
        set(gca, 'XTick', []); % Optionally remove x-axis ticks if desired
        set(gca, 'YTick', []); % Optionally remove y-axis ticks if desired
        saveas(gcf,strcat('../subfigures/Fig4_02132024_ground_truth_exp_',num2str(i_ex)),'epsc');
        close()
        
        
    end
end

%%
if 1
    field_names = fieldnames(methods_struct);
    
    
    for i_ex = 1:length(example_ind_reordered)
        for i_field = 1: length(field_names)
            % methods_struct are older version ordered
            clusters_to_comp = methods_struct.(field_names{i_field}).cluster_array_all{example_ind(i_ex)};
            
            
            b = clusters_to_comp;

            % non_responders are newer version ordered
            non_responder_cluster_drug_regime_codon = non_responders_array{example_ind_reordered(i_ex)};% is this defined by coodn?

            % define the non-responder cluster for each method as the one overlapping with
            % the non-responder cluster the most:
            max_overlap = 0;
            most_overlapping_element = [];

            for i = 1:length(b)
                current_element = b{i};

                % Ensure it's a row vector for comparison
                if isscalar(current_element)
                    current_element = [current_element];
                else
                    current_element = current_element(:)'; % Convert to row vector
                end

                % Find intersection
                overlap = intersect(current_element, non_responder_cluster_drug_regime_codon);

                % Check if this is the maximum overlap
                if length(overlap) > max_overlap
                    max_overlap = length(overlap);
                    most_overlapping_element = current_element;
                end
            end

            non_responder_cluster_elements = most_overlapping_element;% non_responder.(field_names{i_field}){i_ex};
            
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
            figure(1)
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
            
            a(triu(true(size(a)), 1)) = NaN;
            imagesc(a);
            colormap(mycolor_map); % Apply your color map
            caxis([0,1]);
            set(gca, 'box', 'off'); % Turn off the box surrounding the plot
            set(gca, 'XAxisLocation', 'bottom'); % Keep the x-axis (bottom)
            set(gca, 'YAxisLocation', 'left'); % Keep the y-axis (left)
            set(gca, 'XTick', []); % Optionally remove x-axis ticks if desired
            set(gca, 'YTick', []); % Optionally remove y-axis ticks if desired
            saveas(gcf,strcat('../subfigures/Fig4_02132024',field_names{i_field},'_exp_',num2str(i_ex)),'epsc');
            close()
            
        end
        
    end
end
