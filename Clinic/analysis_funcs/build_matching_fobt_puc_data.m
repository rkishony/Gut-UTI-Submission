function cond = build_matching_fobt_puc_data(DrugMapping, cond_row, ...
    POP, Cohorts, PUC, MHS_ResGroups, I_as_RS, DayEdges, RF, RF_Entities, PHEN, PHEN_Entities)
% Create a struct that link FOBT values and PUC values
% PUC are returned as v123count which count, for each sample, drug and time bin, 
% the number of PUCs that are (1) Null (no measurement); (2) Sen; (3) Res

% Sizes:
% nU - Num PUCs
% nD - Num Drugs
% nF - Num FOBT
% nT - Num of time bins
%
% *** Data flow for arranging PUC per FOBT sample pe time ***
%
%  PUC.vnum [nU,all-drugs]  nan: Null, 0: S, 1: I, 2-3: R.
%                           |
%                           | combine_resistance_measurements
%                           | params: MHS_ResGroups
%                           V
%      vnum [nU,nD]         nan: Null, 0: S, 1: I, 2-3: R.
%                           |
%                           | convert_vnum_to_vnan01
%                           | params: I_as_RS
%                           V
%    vnan01 [nU,nD]         nan: Null, 0: Sen, 1: Res
%                           |
%                           | pivot_mhs_res_by_ids_and_time
%                           | params: DayEdges
%                           V
% v123count [nF,nD,nT,3]    counts the number of Null, Sen, Res in each time bin

cond = table2struct(cond_row);
label = cond.type + '_' + cond.tag;

% Choose Resfinder (RF) or Phenotyping (PHEN):
if cond.type == "RF"
    % Resfinder data
    FDATA = RF;
    FDATA_Entities = RF_Entities;
    jFDATA = POP.jRF;
else
    % Phenotypic data
    FDATA = PHEN;
    FDATA_Entities = PHEN_Entities;
    jFDATA = POP.jPHEN;
end

% Restict FOBT population:
cohort = Cohorts.(cond.pop_mask);
jFDATA = jFDATA(cohort);
FDATA = FDATA(jFDATA,:);

% Restrict PUC:
cPUC = PUC(cond.puc_mask(PUC),:);

% Loop over drugs
jDrugs = find(DrugMapping.type == cond.type);
nDrugs = numel(jDrugs);
fobt_val_m = nan(height(FDATA), nDrugs);
nDayBinbs = numel(DayEdges) - 1;
v123count_m = nan(height(FDATA), nDrugs, nDayBinbs, 3);
for ijDrugs = 1:nDrugs
    jDrug = jDrugs(ijDrugs);
    dmap = table2struct(DrugMapping(jDrug,:));

    % For each FOBT, get the pathobiome value (e.g. covergage, num mutations, num colonies):
    fobt_val = get_FOBT_values(FDATA, FDATA_Entities, dmap.entity, dmap.field);
    fobt_val_m(:,ijDrugs) = fobt_val;

    % For each FOBT, get the PUC counts across time:
    vnum_chosen_ab = combine_resistance_measurements(cPUC.vnum(:,MHS_ResGroups.(dmap.ab_group)), 'max', 2);
    vnan01 = convert_vnum_to_vnan01(vnum_chosen_ab, I_as_RS);

    v123count = pivot_mhs_res_by_ids_and_time(...
        FDATA{:,'FOBT_ID'}, cPUC{:,{'FOBT_ID','DateDiff'}}, vnan01, DayEdges);
    v123count_m(:,ijDrugs,:,:) = v123count;

end

cond.label = label;
cond.fobt_val = fobt_val_m;
cond.jFDATA = jFDATA;
cond.v123count = v123count_m;
cond.jDrugs = jDrugs;
cond.nDrugs = numel(jDrugs);

end