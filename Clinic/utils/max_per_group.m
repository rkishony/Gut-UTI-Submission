function grpv = max_per_group(v, labels, ulabels)
if nargin<3
    ulabels = unique(labels);
end
[~, julabels] = ismember(labels, ulabels);
grpv = nan(size(v,1), numel(ulabels));
for i = 1:numel(ulabels)
    m = julabels == i;
    if any(m) 
        grpv(:,i) = max(v(:,m), [], 2);
    end
end

end
