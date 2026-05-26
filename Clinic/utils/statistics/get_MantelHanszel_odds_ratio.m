function OR = get_MantelHanszel_odds_ratio(CMs, dim, remove_dims, pos1, pos2)
% Compute the Mantel-Hanszel Odds Ratio on an array of confusion matrices.
%
%  OR = get_MantelHanszel_odds_ratio(CMs, dim)
%
% where:
%
% CMs [2, 2, n3, n4, ...]    
%
%           an array of [n3, n4, ...] arrays of 2x2 confusion matices.
%
% dim       Dimension(s) accross which to combine the
%           confusion matrices.
%           dim = [] (default): 
%               Do not combine conf matrices.
%               The function just returns the odds ration for each
%               confusion matrix seperately.
%
% Output (example for dim=4):
% OR  is an array of odds-ratio for each confusion matrix pool with size [1, 1, n3, n4, ...], 
% except with singleton at the `dim` positions.
% For example, for dim=4, size OR is [1, 1, n3, 1, ...].
%
%  OR = get_MantelHanszel_odds_ratio(CMs, dim, remove_dims)
% If remove=true, removes singleton dimensions 1, 2 and dim from OR.
%
%  OR = get_MantelHanszel_odds_ratio(CMs, dim, remove_dims, pos1, pos2)
%
% allows using an input array:
%
% CMs [n1, n2, n3, n4, ...]    
%
%           as an array of [n3, n4, ...] arrays of n1xn2 confusion matices
%           (instead of 2x2), with pos1, pos2 as size-2 vectors of indices 
%           specifying the slice in positions 1, 2 for creating a 2x2 array
%           CMs(pos1, pos2, n3, n4, ...)
%

if nargin<2
    dim = [];
end
if nargin<3
    remove_dims = false;
end
if nargin<4 || isempty(pos1)
    assert(size(CMs,1) == 2)
    pos1 = [1 2];
end
if nargin<5 || isempty(pos2)
    assert(size(CMs,2) == 2)
    pos2 = [1 2];
end

% Extract the arrays for the 4 elements of the confusion matrix:
%
%  a  c
%  b  d

cm_cell = num2cell(CMs, 3:ndims(CMs)); % [2, 2] cell aray
[a, b, c, d] = cm_cell{pos1, pos2};   % a, b, c, d are each [1, 1, n3, n4, ...]

% Calculate Mantel-Hanszel (see Wikipedia):
n = a + b + c + d;

ad = a .* d;
bc = b .* c;

if isempty(dim)
    % Regular odds ratio
    OR = ad ./ bc;
else
    % Mantel-Hanszel odds ratio
    M = sum(ad./n, dim);
    D = sum(bc./n, dim);
    OR = M./D;
end

% Remove dims 1 and 2 (the confusion maatrix dims) and the summation dims specified by `dim`
if remove_dims
    rdims = [1 2 dim];
    assert( all(size(OR, rdims) == 1) )
    OR = permute(OR, [setdiff(1:ndims(OR), rdims), rdims]);
end

end
