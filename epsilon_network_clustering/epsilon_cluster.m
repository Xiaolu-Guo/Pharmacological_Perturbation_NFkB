
function cluster_array = epsilon_cluster(pw_dist, epsilon, max_clusters)

    cluster_options = [max_clusters:-1:1];
    ind = 1;

    while ~isempty(cluster_options)
        seed_col = cluster_options(end);
    
        cluster_options(cluster_options == seed_col) = [];
    
        seed = pw_dist(:, seed_col);

        seed_nhbrs = find(seed <= epsilon); %all points that are within epsilon distance of seed (always includes itself)
    
        cluster_array{ind} = seed_nhbrs;
        next_seeds = cluster_array{ind};
        next_seeds(next_seeds == ind) = []; %ind is already part of cluster
        while ~isempty(next_seeds)
        
            seed_col = next_seeds(end);
            seed = pw_dist(:, seed_col);
            seed_nhbrs = find(seed <= epsilon);
            next_seeds(end) = []; %remove last seed from search
        
            for ii = 1:length(seed_nhbrs)
                
                cluster_options(cluster_options == seed_col) = [];
                
                if any(ismember(cluster_array{ind}, seed_nhbrs(ii)))
                    continue
                else               
                    next_seeds = [[next_seeds]; [seed_nhbrs(ii)]];
                    cluster_array{ind} = [cluster_array{ind}; seed_nhbrs(ii)];
                end
            
            end

        end
        ind = ind + 1;
    end
end