function permuted_ids = permute_fobt_ids(ids)
uni_ids = unique(ids);
n = length(uni_ids);
p = randperm(n);
permuted_ids = nan(size(ids));
for i = 1:length(uni_ids)
    if ~isnan(uni_ids(i))
        permuted_ids(ids==uni_ids(i)) = uni_ids(p(i));
    end
end
end