function [CM, uA, uB] = confusion_expand(A, B, uA, uB)
% CONFUSION_EXPAND   Per-slice confusion matrices between A and B over
% sample dim 1
%   [CM, uA, uB] = CONFUSION_EXPAND(A,B) returns an array of size
%     [numel(UA), numel(UB), D2, D3, ...]
%   where UA = unique(A(:)), UB = unique(B(:)), and [D2,D3,...] is the
%   common singleton-expanded size of A and B beyond dims 1-2.
%
%   CM = CONFUSION_EXPAND(A,B,uA,uB) allows specifying the unique values of
%   A and B.
%

% 1) Broadcast to common size
[A,B] = singleton_expand(A,B);
sz = size(A);

% 2) Flatten sample dim vs tasks dims
n      = sz(1);
tasksz = sz(2:end);
P      = prod(tasksz);

% 3) Unique labels + indices
if nargin<3
    uA = unique(A(:));
end
if nargin<4
    uB = unique(B(:));
end
[~, iA] = ismember(A(:), uA);
[~, iB] = ismember(B(:), uB);

nuA = numel(uA);
nuB = numel(uB);

% 4) Build group indices for each slice and accumulate
grp  = repmat(1:P, n, 1);
CMf  = accumarray([iA, iB, grp(:)], 1, [nuA, nuB, P]);

% 5) Reshape to [uA x uB x tasksz]
CM = reshape(CMf, [nuA, nuB, tasksz]);
end
