function OR = get_MantelHanszel_odds_ratio_for_F3xV4(CMs, dim, remove_dims, pos1, pos2)
% Calls get_MantelHanszel_odds_ratio with default param set for the confusion matrix
% array of FOBT SEN/RES association with PUC V4 [NONE, NULL, SEN, RES]
% array, and assuming drugs are in dimension 3:
%
% CMs  [2, 4, nDrugs, ...]
%

if nargin<2
    dim = 3; % Sum over the 'Drug' dimension.
end
if nargin<3
    remove_dims = true;
end
if nargin<4
    pos1 = 2:3;  % SEN/RES of FOBT
end
if nargin<5
    pos2 = 3:4;  % SEN/RES of PUC
end

OR = get_MantelHanszel_odds_ratio(CMs, dim, remove_dims, pos1, pos2);

end