load("non_responders_inds_011924.mat") %Old non-responder threshold

%Store non-responder indices in their corresponding drug regime
%cluster_array cell
i_start = 1;
i_end = 15;
num_conditions = 15;
drug_regime = 1;
while i_end <= 3000
    non_responders_array{drug_regime} = [];
    for i = 1:length(no_response_row)
        if any(i_start:i_end == no_response_row(i))
    
            non_responders_array{drug_regime} = ...
            [non_responders_array{drug_regime}; no_response_row(i)];
        end
    end
    
    i_start = i_start + 15;
    i_end = i_end + 15;
    drug_regime = drug_regime + 1;
    
end

%% 
%Converting 1:3000 indices to 1:15 indices within each cluster_array cell
%based on its stimulus condition
i_start = 1;
i_end = 15;

for drug_regime = 1:length(non_responders_array)
    
    non_responders_new_ind_array{drug_regime} = [];
    ind_range(1, :) = [i_start:i_end];
    ind_range(2, :) = [1:15];
    if ~isempty(non_responders_array{drug_regime})
        for ind = 1:length(non_responders_array{drug_regime})
            
            old_ind = non_responders_array{drug_regime}(ind);
            new_ind = ind_range(2, ind_range(1, :) == old_ind);
            non_responders_new_ind_array{drug_regime} = ...
                [non_responders_new_ind_array{drug_regime}, new_ind];
        end
        
        
    end
    
    i_start = i_start + 15;
    i_end = i_end + 15;
end

%%
non_responders_array = non_responders_new_ind_array;

save("non_responders_per_drug_regime.mat", "non_responders_array")