function [POP, FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN] = ...
    merge_FOBT_with_MHS(FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN, Timespan)

% Add demographics to FOBT:
FOBT = addMatch(FOBT, IDS, 'RandomID', {'Gender','age'});

% Closest PUC after each FOBT
FOBT = addClosest(FOBT, SMP);

% Find if we have any PUC (0-365):
is_following_fobt_bac = @(fobt,bac) inrange(bac.SampleDate - fobt.FOBT_Date, Timespan.BAC_FOBT_inclusion);
FOBT.anyPUC = countMatches(FOBT, BAC, 'FOBT_ID', is_following_fobt_bac, true);
if Timespan.BAC_FOBT_inclusion == 0
    assert(all(FOBT.anyPUC == (FOBT.nearest_SMP_j(:,2)<=Timespan.BAC_FOBT_inclusion(2))))
end

% Find coinciding FOBT and BAC:
is_coinciding_fobt_bac = @(fobt,bac) inrange(bac.SampleDate - fobt.FOBT_Date, Timespan.BAC_FOBT_coinciding);
[FOBT.CoincidingPUC, BAC.CoincidingFOBT] = countMatches(FOBT, BAC, 'FOBT_ID', is_coinciding_fobt_bac, true);

% POP is defined as a subset of FOBT with unique patients and followup PUCs:
POP = FOBT(FOBT.anyPUC & FOBT.FOBTisFirst,:);

% Identify for each PUC whether it has the same bacterial species as the FOBT-coinciding PUC:
is_bac_same_species_as_coinciding_bac = @(bac1, bac2) bac2.CoincidingFOBT & bac1.Bacteria == bac2.Bacteria;
BAC.isSameSpeciesAsCoinciding = countMatches(BAC, BAC, 'FOBT_ID', is_bac_same_species_as_coinciding_bac, true);
%%% TODO: Can be statistically more powerful to set the SMP to the BAC that
%         is not matching (instead of discarding).
SMP.isSameSpeciesAsCoinciding = countMatches(SMP, BAC(BAC.isSameSpeciesAsCoinciding,:), 'UTI_Drisha', [], true);

% Add demographics to RESFDR and link with FOBT and POP
[RESFDR, ~, POP.jRF] = merge_with_FOBT_and_POP(RESFDR, FOBT, POP, 'RESFDR');

% Add demographics to PHEN and link with FOBT and POP
[PHEN, ~, POP.jPHEN] = merge_with_FOBT_and_POP(PHEN, FOBT, POP, 'PHEN');
j = POP.jPHEN; j = j(j>0);
assert(all(PHEN.jFOBT(j)>0 & ~PHEN.isrep(j) & ~PHEN.is_neg_ctrl(j)))

% Add UTI diag to BAC and SMP:
is_matching_bac_diag = @(bac,diag) inrange(diag.date_diagnosis - bac.SampleDate, Timespan.DIAG_BAC_coinciding);
BAC.isuti = countMatches(BAC, DIAG(DIAG.UTI,:), 'FOBT_ID', is_matching_bac_diag) > 0;
SMP = addMatch(SMP, BAC, {'UTI_Drisha', 'FOBT_ID'}, {'isuti', 'CoincidingFOBT'});

end


function [TBL, FOBT_jTBL, POP_jTBL] = merge_with_FOBT_and_POP(TBL, FOBT, POP, TBL_name)
TBL = addMatch(TBL, POP, 'FOBT_ID', {'FOBT_Date','RandomID','FOBTisFirst','CoincidingPUC', 'Gender', 'age'});

[~,TBL.jFOBT] = ismember(TBL.FOBT_ID, FOBT.FOBT_ID);
missing_in_FOBT = TBL.jFOBT==0 & TBL.isFOBT;
if any(missing_in_FOBT)
    fprintf('WARNING: Removing from %s %d FOBT_ID`s not found in FOBT table\n', TBL_name, sum(missing_in_FOBT))
    TBL(missing_in_FOBT, :) = [];
end
assert(all((TBL.jFOBT>0) == TBL.isFOBT))

[~,TBL.jPOP] = ismember(TBL.FOBT_ID, POP.FOBT_ID);
assert(all((TBL.jPOP>0) == (TBL.isFOBT & TBL.FOBTisFirst)))

[~,POP_jTBL] = ismember(POP.FOBT_ID, TBL.FOBT_ID);
assert(sum(TBL.isrep(POP_jTBL(POP_jTBL>0)))==0)

[~,FOBT_jTBL] = ismember(POP.FOBT_ID, TBL.FOBT_ID);

TBL.isFOBT = [];

end


function FOBT = addClosest(FOBT, SMP)
h = height(FOBT);
FOBT.nearest_SMP_j = zeros(h, 2);
FOBT.nearest_SMP_t = zeros(h, 2);

for jF = 1:height(FOBT)
    jS = find(SMP.FOBT_ID == FOBT.FOBT_ID(jF));
    assert(all(SMP.DateDiff(jS) == sort(SMP.DateDiff(jS))))
    is_pos = SMP.DateDiff(jS) >= 0;
    f_pos = find(is_pos, 1);
    if isempty(f_pos)
        SMP_j2 = 0;
        SMP_t2 = inf;
    else
        SMP_j2 = jS(f_pos);
        SMP_t2 = SMP.DateDiff(SMP_j2);
    end

    f_neg = find(~is_pos, 1, "last");
    if isempty(f_neg)
        SMP_j1 = 0;
        SMP_t1 = -inf;
    else
        SMP_j1 = jS(f_neg);
        SMP_t1 = SMP.DateDiff(SMP_j1);
    end
    FOBT.nearest_SMP_j(jF,:) = [SMP_j1, SMP_j2];
    FOBT.nearest_SMP_t(jF,:) = [SMP_t1, SMP_t2];
end


end