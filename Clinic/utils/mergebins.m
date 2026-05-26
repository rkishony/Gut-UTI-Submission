function out = mergebins(count, dim, edges, new_edges)
%MERGEBINS Merge histogram/bin counts into coarser bins.
%   OUT = MERGEBINS(COUNT, DIM, EDGES, NEW_EDGES) merges counts in COUNT
%   along dimension DIM from EDGES to NEW_EDGES. Assumes NEW_EDGES is a
%   subset of EDGES and bin boundaries coincide (no partial overlaps).

edges = edges(:); new_edges = new_edges(:);
assert(isscalar(dim) && dim>=1 && dim<=ndims(count), 'dim out of range');
assert(all(ismember(new_edges, edges)), 'NEW_EDGES must be subset of EDGES');

nOld = numel(edges) - 1;
sz = size(count);
assert(sz(dim) == nOld, 'Size along DIM must equal numel(EDGES)-1');

% Move DIM to front and flatten remainder
p = [dim, setdiff(1:ndims(count), dim)];
cntP  = permute(count, p);
tail  = size(cntP); tail = tail(2:end);
cnt2D = reshape(cntP, nOld, []);

% Group ids for old bins using boundary markers from ismember
mark = ismember(edges, new_edges);
assert(mark(1) && mark(end), 'NEW_EDGES must span EDGES');
g    = cumsum(mark(1:end-1));
nNew = sum(mark) - 1;

% Accumulate with 2-D subs (works for multi-col cnt2D)
m    = size(cnt2D, 2);
subs = [repmat(g, m, 1), repelem((1:m).', nOld, 1)];
new2D = accumarray(subs, cnt2D(:), [nNew, m], @sum, 0);

% Restore shape and dimension order
outP = reshape(new2D, [nNew, tail]);
out  = ipermute(outP, p);
end
