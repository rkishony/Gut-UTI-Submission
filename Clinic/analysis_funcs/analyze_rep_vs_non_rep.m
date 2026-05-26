function [tot_num_replicates, tot_unique_ids_with_replicates, dists, p_vals] = analyze_rep_vs_non_rep(data, ids)
data = double(data);
n = size(data,1);

[~,i] = unique(ids);
is_chosen = false(n,1);
is_chosen(i) = true;
is_chosen2 = is_chosen & is_chosen';

is_rep = count_occurrences(ids) > 1;
is_rep2 = is_rep & is_rep';

tot_num_replicates = sum(is_rep);
tot_unique_ids_with_replicates = numel(unique(ids(is_rep)));

warning('off')
dist = pdist2(data, data, 'hamming');
warning('on')

same_id = ids == ids';
valid = triu(ones(n),1);

masks{1} =  same_id & valid; % replicates
masks{2} = ~same_id & is_rep2 & is_chosen2 & valid; % different ids within the rep ids  (one sample per patient)
masks{3} = ~same_id & is_chosen2 & valid; % all different id's (one sample per patient)

dists = cellfun(@(m)dist(m)', masks, 'UniformOutput',false);

[p_vals.p12,~,~] = ranksum(dists{1}, dists{2});
[p_vals.p13,~,~] = ranksum(dists{1}, dists{3});
[p_vals.p23,~,~] = ranksum(dists{2}, dists{3});

end