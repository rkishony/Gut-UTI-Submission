function [CM, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile, fit_nominal, fit_bs_prctile] = ...
    calc_F3_V4_association_across_time_and_drugs( ...
    F3, V4, Days, prctiles, randomization, num_perm, num_bootstrap)

if nargin<6 || isempty(num_perm)
    num_perm = 10000;
end
if nargin<7 || isempty(num_bootstrap)
    num_bootstrap = 10000;
end

% Calculate confusion matrices (CM) between F3 (FOBT) and V4 (PUC):

[CM, CM_perm, CM_bs] = calc_confmat_over_drugs_and_time_with_permutations(...
    F3, V4, randomization, num_perm, num_bootstrap);
% CM      : [3, 4, nD, nT     ]
% CM_perm : [3, 4, nD, nT, nP ]
% CM_bs   : [3, 4, nD, nT, nBS]

% Calculate cross-drugs Mantel-Hanszel odds ratios:
OR = get_MantelHanszel_odds_ratio_for_F3xV4(CM);
OR_perm = get_MantelHanszel_odds_ratio_for_F3xV4(CM_perm);
OR_bs = get_MantelHanszel_odds_ratio_for_F3xV4(CM_bs);

% Calculate P-value based on permutation analysis:
OR_perm_P_val = mean(OR_perm >= OR, 2);

% Calculate Confidence Interval IC90:
OR_perm_prctile = prctile(OR_perm, prctiles, 2);
OR_bs_prctile = prctile(OR_bs, prctiles, 2);

% Calculate linear fit of log(OR) vs time:

log_OR_bs = log(OR_bs);
log_OR_bs(log_OR_bs==-inf) = 1e-3;
log_OR_bs(log_OR_bs==+inf) = 1e+3;

fit_bs = zeros(num_bootstrap, 2);
for i = 1:num_bootstrap
   fit_bs(i,:) = polyfit(Days, log_OR_bs(:,i), 1);
end
fit_bs_prctile = prctile(fit_bs, prctiles, 1);

fit_nominal = polyfit(Days, log(OR), 2);

end