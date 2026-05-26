function v123count = pivot_mhs_res_by_ids_and_time(unique_ids, ids_days, vnan01, day_edges)
% Convert a list of PUC events into an array of events per timebin
%
% unique_ids    [nF, 1]         Vector of unique IDs
% ids_days      [nU, 2]         Array of PUC events, indicating patient ID and
%                               time (days since FOBT)
% vnan01        [nU, 1]         Vector of PUC with nan fe Null measurement, 
%                               0 for Sen and 1 for Res
% day_edges     [1, nT]         Vector indicating the edges of the time bins
%                               into which PUC are counted.
%
% Returns
%
% v123count     [nF, nT-1, 3]   Array of counts of PUCs per time bin, 
%                               1: Null, 2: Sen, 3: Res.

assert(iscolumn(unique_ids))
assert(iscolumn(vnan01))

nIDs = numel(unique_ids);
nBinbs = numel(day_edges) - 1;
assert(size(ids_days,1) == size(vnan01,1))

% Convert [nan,0,1] -> [1,2,3]
v123 = vnan01;
v123(isnan(v123)) = -1;
v123 = v123 + 2;

jTime = discretize(ids_days(:,2), day_edges);
% assert(~any(isnan(jTime)));

[~,jID] = ismember(ids_days(:,1), unique_ids);

subs = [jID, jTime, v123];
ok = all(subs>0, 2);
if any(~ok)
    % fprintf('Removing %d non-matching measurements\n', sum(~ok))
end

v123count = accumarray(subs(ok,:), 1, [nIDs, nBinbs, 3]);

end