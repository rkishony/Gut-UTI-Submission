function sortIdx = sort_rows_by_cluster(data, metric, linkageMethod)
if nargin < 2 || isempty(metric)
    metric = 'euclidean';
end
if nargin < 3 || isempty(linkageMethod)
    linkageMethod = 'average';
end

Y = pdist(data, metric);
Z = linkage(Y, linkageMethod);
sortIdx = optimalleaforder(Z, Y);
end
