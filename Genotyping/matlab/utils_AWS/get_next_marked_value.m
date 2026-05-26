function [value, idx] = get_next_marked_value(A, M, dim)
% A     : input array
% M     : logical array same size as A (true = aamarkedaa)
% dim   : dimension to search along
%
% idx   : same size as A, giving the 1-based index
%         along dim of the next marked entry (0 if none)
% value : same size as A, giving A(idx) (NaN if idx==0)

if isempty(M)
    M = A>0;
elseif isscalar(M) && ~islogical(M)
    M = A>M;
end

assert(all(size(M)==size(A)))
assert(islogical(M))

% reshape so dim is last
sz   = size(A);
n    = sz(dim);
m    = numel(A)/n;
perm = [setdiff(1:ndims(A),dim), dim];

Ap   = permute(A, perm);
A2   = reshape(Ap, [m, n]);

Mp   = permute(M, perm);
M2   = reshape(Mp, [m, n]);

% compute next-index along rows
idxs   = 1:n;
revIdx = M2(:, end:-1:1) .* idxs;
revC   = cummax(revIdx, 2);
posRev = revC(:, end:-1:1);
idx2   = (posRev~=0) .* (n+1 - posRev);  % 0 where no next

% gather values
value2       = nan(size(A2));
rows         = (1:m)';
linIdxMatrix = (idx2 - 1) * m + rows(:,ones(1,n));
valid        = idx2 > 0;
value2(valid)= A2(linIdxMatrix(valid));

% reshape & undo permutation
idx   = ipermute( reshape(idx2,   size(Ap)), perm );
value = ipermute( reshape(value2, size(Ap)), perm );
end
