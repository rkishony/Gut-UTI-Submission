function [CM, CM_perm, CM_bs] = calc_confmat_over_drugs_and_time_with_permutations(...
    F3, V4, randomization, num_perm, num_bs)
% Calculate nominal, permuted and boostrapped confusion matrices 
% between FOBT and PUC over drugs and time.
%
% CM = calc_confmat_over_drugs_and_time_with_permutations(F3, V4)
%
% F3          [nF, nD]        Array of 1-3 for nan/sensitive/resistant FOBT for
%                             each patient and each drug.
%
% V4          [nF, nD, nT]    Array of resistance type 1-4 (NONE, NULL, SEN, RES)
%
%
% Returns:
%
% CM        [3, 4, nD, nT]  
%
%  CM provides for each drug and each timebin, a 3x4 confusion matrix:
%
%                    ------- PUC -------
%          FOBT      NONE NULL SEN  RES
%          nan (1)
%          sen (2)   -    -    -    -
%          res (3)   -    -    -    -
% 
% The columns indicate:
%    1 - NONE No PUCs.
%    2 - NULL Only PUCs with no measurements.
%    3 - SEN  At least 1 sensitive PUC and no resistant PUCs.
%    4 - RES  At least 1 resistant PUC.
%
%
%  [CM, CM_perm] = calc_confmat_over_drugs_and_time_with_permutations(___, randomization, perm_num)
% Also permutes the rows of F3 `num_perm` times and returns:
%
% CM_perm   [3, 4, nD, nT, num_perm]
%
%  randomization:   'complete' - permute all the set
%                   'restricted' - permute only within the measured set (namely, F3=2,3 and V4=3,4)
%
%  [CM, CM_perm, CM_bs] = calc_confmat_over_drugs_and_time_with_permutations(___, perm_num, bs_num)
% Also bootsrap the samples `n_bs` times and returns:
%
% CM_bs     [3, 4, nD, nT, num_bs]
% 

if nargin<4 || isempty(num_perm)
    num_perm = 0;
end
if nargin<5 || isempty(num_bs)
    num_bs = 0;
end

% Calculate nominal:
CM = confusion_expand_f3v4(F3, V4);  % [2, 4, nD, nT]

n = size(F3,1);

switch randomization
    case 'complete'
        inds_to_randomize = 1:n;
    case 'restricted'
        inds_to_randomize = find(all(F3>=2,2));
    otherwise
        error('Unknown randomization, %s', randomization)
end
n_rand = numel(inds_to_randomize);

% Calculate permutations:
if nargout>1
    fprintf('Permutation (%d) ', num_perm)
    CM_perm = nan([size(CM), num_perm]);
    for i = 1:num_perm
        if ~mod(i, round(num_perm) / 20)
            fprintf('.')
        end
        ind = randperm(n_rand);
        F3_perm = F3;
        F3_perm(inds_to_randomize,:) = F3(inds_to_randomize(ind),:);
        CM_perm(:,:,:,:,i) = confusion_expand_f3v4(F3_perm, V4);
    end
    fprintf('\n')
end

% Calculate bootstraps:
if nargout>2
    fprintf('Bootstrap   (%d) ', num_bs)
    CM_bs = nan([size(CM), num_bs]);
    for i = 1:num_bs
        if ~mod(i, round(num_bs) / 20)
            fprintf('.')
        end
        ind = randi(n_rand, 1, n_rand);
        F3_bs = F3;
        F3_bs(inds_to_randomize,:) = F3(inds_to_randomize(ind),:);
        V4_bs = V4;
        V4_bs(inds_to_randomize,:,:) = V4(inds_to_randomize(ind),:,:);
        CM_bs(:,:,:,:,i) = confusion_expand_f3v4(F3_bs, V4_bs);
    end
    fprintf('\n')
end

end



function cm = confusion_expand_f3v4(f3, v4)
cm = confusion_expand(f3, v4, 1:3, 1:4);
end
