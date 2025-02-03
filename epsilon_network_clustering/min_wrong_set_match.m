function [wrong_num_min,set_matches_min] = min_wrong_set_match(set1,set2)

set_match = clustering_compare(set1,set2);

for i_match = 1:length(set_match)
    wrong_num(i_match) = set_match{i_match}.set1_wrong_match_ele_num;
end

[wrong_num_min,wrong_num_min_id] = find(wrong_num == min(wrong_num));

wrong_num_min = wrong_num_min(1);
wrong_num_min_id = wrong_num_min_id(1);

set_matches_min = set_match{wrong_num_min_id}.set12_index_match;
wrong_num_min = set_match{wrong_num_min_id}.set1_wrong_match_ele_num;
end