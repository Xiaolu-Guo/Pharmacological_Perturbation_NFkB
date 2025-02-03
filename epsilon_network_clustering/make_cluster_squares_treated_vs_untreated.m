function make_cluster_squares_treated_vs_untreated(cluster_array, non_responders_array, num_row, num_col, subplot_order, figure_path)
figure('units','inch','position',[0,0,3.045,3.045])
my_color_map = [0,0,0;
    1,0,0;
    1,1,1];
%figure()
cc = 1;

for i = 1:length(cluster_array)
    
    cluster_square = ones(15, 15);
    for m = 1:15
        ut_col = cluster_array{i}{m};
        if numel(ut_col) == 0
            cluster_square(:, m) = 1;
        else
            for t_ind = 1:length(ut_col)
                n = cluster_array{i}{m}(t_ind);
                cluster_square(n, m) = 0;
                if m == 4 %non-responder LPS LSD
                    cluster_square(n, m) = 0.5;
                else
                    cluster_square(n, m) = 0;
                end
            end
        end
    end

    subplot(num_row, num_col, subplot_order(cc))

    imagesc(cluster_square)
    caxis([0, 1])
    %'TickLength', [0 0],
    set(gca,'TickLength', [0 0], 'ytick', [1:15], 'yticklabel', [],...
    'xtick', [1:15], 'xticklabel', [], 'FontName', 'Arial',...
    'FontSize', 11)
    %set(gcf, 'Position', [100 100 800 800])
    %set(gcf, 'PaperUnits', 'inches');
    %set(gcf, 'PaperSize', [7 7]);
    xtickangle(270)
    colormap(my_color_map)
    
    ax=gca;
    ax.LineWidth = 1;
    clear cluster_square
    cc = cc + 1;
    
    
    saveas(gcf, figure_path, 'png')
end
end