function [b_same, b_bs, nprofs, nperm, by_chance_in_nperm] = phenotypic_distance_analysis(IsolateResDataTable, Abs, nperm)

if nargin<3 || isempty(nperm)
    nperm = 1e4;
end

%% unpack table into vectors

isMeasured = is_measured(IsolateResDataTable);

valid = is_valid(IsolateResDataTable,isMeasured);

IsolateResDataTable = IsolateResDataTable(valid,:);
resistance_maps = IsolateResDataTable{:, Abs};
fobt_id = IsolateResDataTable.FOBT_ID;
uti_id = IsolateResDataTable.UTI_ID;
is_uti = IsolateResDataTable.UTI_ID>0 & ~IsolateResDataTable.is_neg_ctrl;
is_fobt = ~is_uti & fobt_id>0 & ~isnan(fobt_id) & ~IsolateResDataTable.is_neg_ctrl;

% [~,j_experiment_code] = find(table2array(IsolateResDataTable(:,{'IsolateResA','IsolateResB','IsolateResC'})));

%% unique lists:
% unique patients in dataset
fpatients = unique(fobt_id(is_fobt));
assert(~any(fpatients==0))

% unique utis in dataset
utis_measured = uti_id(is_uti>0);
uti_patients = fobt_id(is_uti>0);
UFpatients = unique(intersect(uti_patients, fpatients));
uti_w_fobt = utis_measured(ismember(uti_patients,UFpatients));
uti_w_fobt_patients = uti_patients(ismember(uti_patients, UFpatients));
[uti_ids, ia, ~] = unique(uti_w_fobt);
uti_patient_ids = uti_w_fobt_patients(ia);

% how many unique fobt profiles there seem to be:
nprofs = nan(size(fpatients));
personal_resmaps = cell(size(fpatients));
for f = 1:length(fpatients)
    fi = fpatients(f);
    ind = find(fobt_id==fi & is_fobt);
    personal_resmaps{f} = resistance_maps(ind,:);
    nprofs(f) = size(unique(personal_resmaps{f},'rows'),1);
end



%%
% uti to same-patient vs any-other pathobiome distances:
compute_dist2other_too = true;
permute_patients = false;

[mindists_same, mindists_other] = computeUTI2SameAndOther(resistance_maps, fobt_id, uti_id, ...
    is_fobt, is_uti, uti_ids, uti_patient_ids, compute_dist2other_too);

% uti to "same-patient" pathobiome -- but permuted:
compute_dist2other_too = false;
permute_patients = true;
perms_mindists_same = nan(length(mindists_same),nperm);
%%
for perm = 1:nperm
    permuted_patient_ids = permute_fobt_ids(uti_patient_ids);
    perms_mindists_same(:,perm) = computeUTI2SameAndOther(resistance_maps, fobt_id, uti_id,...
        is_fobt, is_uti, uti_ids, permuted_patient_ids, compute_dist2other_too);
end



%%

binning = -0.5:length(Abs)+0.5;
hc_same  = histcounts(mindists_same,  binning);
hc_other = histcounts(mindists_other, binning);
b_same  = hc_same /sum(hc_same);
b_other = hc_other/sum(hc_other);

for perm = 1:nperm
    hc_bs(:, perm) = histcounts(perms_mindists_same(:,perm),binning);
    b_bs(:, perm) = hc_bs(:,perm)/sum(hc_bs(:,perm));
end
by_chance_in_nperm = sum(b_bs(1,:)>=b_same(1));
end

