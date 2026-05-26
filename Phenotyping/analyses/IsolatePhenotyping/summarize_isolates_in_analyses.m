function A = summarize_isolates_in_analyses(T)
% BUILD_ANALYSIS_STRUCT (Updated Jul-2025)
%
% Rules:
%    Measured = not FromLawnOrEmpt
%    Valid    = Measured AND grew in MC (no_growth_in_mc == false)
%    PhenProfileMatching patients = those with >=1 valid FOBT AND >=1 valid UTI isolate
%    Isolate key = (isolate_id, fobt_id, experiment)

keyfun = @(iso,pid,exp) unique(string(iso) + "_" + string(pid) + "_" + string(exp), 'stable');

% define masks
isMeasured = is_measured(T);
isValid    = is_valid(T, isMeasured);
isFOBT     = T.UTI_ID == 0;
isUTI      = T.UTI_ID > 0;
patient_id = T.FOBT_ID;

% 3. experiment column
expNames = unique(T.expname);

% PhenProfileMatching paired patients
validFOBTpatients = unique(patient_id(isFOBT & isValid));
validUTIpatients = unique(patient_id(isUTI  & isValid));
pairedPatients   = intersect(validFOBTpatients, validUTIpatients);

% for MatchPhenGen
maskB        = T.expname == "IsolateResB";
patientsInB  = unique(patient_id(maskB));

% stash helper
stash = @(rowsFOBT,rowsUTI,expCol) deal( ...
    struct( ...
      'initial_fobt_isolate_ids' , keyfun(rowsFOBT.iso , rowsFOBT.pid , expCol), ...
      'valid_fobt_isolate_ids'   , keyfun(rowsFOBT.isoV, rowsFOBT.pidV, expCol), ...
      'initial_uti_isolate_ids' , keyfun(rowsUTI.iso , rowsUTI.cult , expCol), ...
      'valid_uti_isolate_ids'   , keyfun(rowsUTI.isoV, rowsUTI.cultV, expCol), ...
      'initial_fobt_patient_ids' , unique(rowsFOBT.pid), ...
      'valid_fobt_patient_ids'   , unique(rowsFOBT.pidV), ...
      'initial_uti_patient_ids' , unique(rowsUTI.pid), ...
      'valid_uti_patient_ids'   , unique(rowsUTI.pidV), ...
      'initial_uti_culture_ids' , unique(rowsUTI.cult), ...
      'valid_uti_culture_ids'   , unique(rowsUTI.cultV) ));

for k = 1:numel(expNames)
    thisExp = expNames{k};
    colmask = logical(T.expname == thisExp);

    % base masks
    baseFOBT = colmask & isFOBT & isMeasured;
    baseUTI = colmask & isUTI & isMeasured;

    % store common parts
    rowsFOBT.iso   = T.isolate_id(baseFOBT);
    rowsFOBT.pid   = patient_id(baseFOBT);
    rowsFOBT.isoV  = T.isolate_id(baseFOBT & isValid);
    rowsFOBT.pidV  = patient_id(baseFOBT & isValid);

    rowsUTI.iso   = T.isolate_id(baseUTI);
    rowsUTI.pid   = patient_id(baseUTI);
    rowsUTI.cult  = T.UTI_ID(baseUTI);
    rowsUTI.isoV  = T.isolate_id(baseUTI & isValid);
    rowsUTI.pidV  = patient_id(baseUTI & isValid);
    rowsUTI.cultV = T.UTI_ID(baseUTI & isValid);

    % Isolate Level Phen Diversity
    A.hist_num_phen_res_fobt_iso_profiles_per_pat.exp.(thisExp) = stash(rowsFOBT, struct('iso',[],'pid',[], ...
                                               'cult',[],'isoV',[],'pidV',[],'cultV',[]), thisExp);
    
    % Phen Profile Matching
    maskPhenProfileMatching = ismember(patient_id, pairedPatients);
    rowsFOBTE = filter_rows(rowsFOBT, maskPhenProfileMatching(baseFOBT), maskPhenProfileMatching(baseFOBT & isValid));
    rowsUTIE = filter_rows(rowsUTI, maskPhenProfileMatching(baseUTI), maskPhenProfileMatching(baseUTI & isValid));
    A.hist_phen_profiles_match_same_vs_diff_pat.exp.(thisExp) = stash(rowsFOBTE, rowsUTIE, thisExp);

    % MatchPhenGen
    maskMatchPhenGen = ismember(patient_id, patientsInB);
    rowsFOBTb = filter_rows(rowsFOBT, maskMatchPhenGen(baseFOBT), maskMatchPhenGen(baseFOBT & isValid));
    A.MatchPhenGen.exp.(thisExp) = stash(rowsFOBTb, struct('iso',[],'pid',[], ...
                                                'cult',[],'isoV',[],'pidV',[],'cultV',[]), thisExp);
end

% Add total counts per analysis
fnAnal = fieldnames(A);
for a = 1:numel(fnAnal)
    ANAL = fnAnal{a};
    exps = fieldnames(A.(ANAL).exp);

    initFOBTiso = []; valFOBTiso = [];
    initUTIiso = []; valUTIiso = [];
    initFOBTpat = []; valFOBTpat = [];
    initUTIpat = []; valUTIpat = [];
    initUTIcul = []; valUTIcul = [];

    for e = 1:numel(exps)
        leaf = A.(ANAL).exp.(exps{e});
        initFOBTiso = [initFOBTiso ; leaf.initial_fobt_isolate_ids];
        valFOBTiso  = [valFOBTiso  ; leaf.valid_fobt_isolate_ids];
        initUTIiso = [initUTIiso ; leaf.initial_uti_isolate_ids];
        valUTIiso  = [valUTIiso  ; leaf.valid_uti_isolate_ids];
        initFOBTpat = [initFOBTpat ; leaf.initial_fobt_patient_ids];
        valFOBTpat  = [valFOBTpat  ; leaf.valid_fobt_patient_ids];
        initUTIpat = [initUTIpat ; leaf.initial_uti_patient_ids];
        valUTIpat  = [valUTIpat  ; leaf.valid_uti_patient_ids];
        initUTIcul = [initUTIcul ; leaf.initial_uti_culture_ids];
        valUTIcul  = [valUTIcul ; leaf.valid_uti_culture_ids];
    end

    A.(ANAL).total.initial_fobt_isolates   = numel(unique(initFOBTiso));
    A.(ANAL).total.valid_fobt_isolates     = numel(unique(valFOBTiso));
    A.(ANAL).total.initial_uti_isolates    = numel(unique(initUTIiso));
    A.(ANAL).total.valid_uti_isolates      = numel(unique(valUTIiso));
    A.(ANAL).total.initial_fobt_patients   = numel(unique(initFOBTpat));
    A.(ANAL).total.valid_fobt_patients     = numel(unique(valFOBTpat));
    A.(ANAL).total.initial_uti_patients    = numel(unique(initUTIpat));
    A.(ANAL).total.valid_uti_patients      = numel(unique(valUTIpat));
    A.(ANAL).total.initial_uti_cultures    = numel(unique(initUTIcul));
    A.(ANAL).total.valid_uti_cultures      = numel(unique(valUTIcul));
end
end

function rowsOut = filter_rows(rowsIn, maskInit, maskValid)
rowsOut.iso   = rowsIn.iso(maskInit);
rowsOut.pid   = rowsIn.pid(maskInit);
if isfield(rowsIn, 'cult')
    rowsOut.cult = rowsIn.cult(maskInit);
else
    rowsOut.cult = [];
end
rowsOut.isoV  = rowsIn.isoV(maskValid);
rowsOut.pidV  = rowsIn.pidV(maskValid);
if isfield(rowsIn, 'cultV')
    rowsOut.cultV = rowsIn.cultV(maskValid);
else
    rowsOut.cultV = [];
end
end