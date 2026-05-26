function S = summarize_isolates_in_experiments(T)
% SUMMARIZE_ISOLATES  v6  (Jul-2025)
%
%   Measured = not FromLawnOrEmpt
%   Valid    = Measured AND grew in MC (no_growth_in_mc == false)
%   patient_id = fit_id for every row
%   Isolate uniqueness = (isolate_id , patient_id , experiment)
%   NEW: S now contains grand-total patient counts across all experiments:
%        S.total_fobt_patients_overall
%        S.total_uti_patients_overall
%        S.total_patients_overall      (union of both)

keyfun = @(iso,id) unique(string(iso) + "_" + string(id), 'stable');

res = get_resistance_matrix(T);

% ------------------------------------------------------------
% 1. MEASURED rows only
isMeasured = is_measured(T);

% 2. growth in MC
isValid = is_valid(T, isMeasured);

% 3. experiment column
expNames = unique(T.expname);

% 4. patient_id = fit_id
patient_id = T.FOBT_ID;

% ------------------------------------------------------------
% 5. loop per experiment
S.total_fobt_isolates = 0;
S.total_uti_isolates = 0;
S.valid_fobt_isolates = 0;
S.valid_uti_isolates = 0;

fobtPatientSet = [];   % will hold union of FIT patients
utiPatientSet = [];   % will hold union of UTI patients
fobtPatientSetV = [];
utiPatientSetV = [];
utiCultureSet = [];
utiCultureSetV = [];
for k = 1:numel(expNames)
    thisExp = expNames{k};
    present = logical(T.expname == thisExp);

    isFIT = T.UTI_ID == 0;
    isUTI = T.UTI_ID > 0;

    % ------------ FIT ------------
    maskFIT      = present & isFIT & isMeasured;
    maskFITvalid = maskFIT & isValid;

    keysTotFIT = keyfun(T.isolate_id(maskFIT), patient_id(maskFIT));
    keysValFIT = keyfun(T.isolate_id(maskFITvalid), patient_id(maskFITvalid));

    patTotFIT  = unique(patient_id(maskFIT));
    patValFIT  = unique(patient_id(maskFITvalid));

    % ------------ UTI ------------
    maskUTI      = present & isUTI & isMeasured;
    maskUTIvalid = maskUTI & isValid;

    keysTotUTI = keyfun(T.isolate_id(maskUTI), T.UTI_ID(maskUTI));
    keysValUTI = keyfun(T.isolate_id(maskUTIvalid), T.UTI_ID(maskUTIvalid));

    patTotUTI  = unique(patient_id(maskUTI));
    patValUTI  = unique(patient_id(maskUTIvalid));

    % ------------ UTI culture IDs ------------
    cultTot = unique(T.UTI_ID(maskUTI));
    cultTot = cultTot(~isnan(cultTot) & cultTot~=0); % just to make sure

    cultVal = unique(T.UTI_ID(maskUTIvalid));
    cultVal = cultVal(~isnan(cultVal) & cultVal~=0);

    % ------------ stash per-experiment ------------
    e = struct;
    e.total_fit_isolates  = numel(keysTotFIT);
    e.valid_fit_isolates  = numel(keysValFIT);
    e.total_uti_isolates  = numel(keysTotUTI);
    e.valid_uti_isolates  = numel(keysValUTI);

    e.total_fit_patient_ids = patTotFIT;
    e.valid_fit_patient_ids = patValFIT;
    e.total_uti_patient_ids = patTotUTI;
    e.valid_uti_patient_ids = patValUTI;

    e.total_fit_patients  = numel(patTotFIT);
    e.valid_fit_patients  = numel(patValFIT);
    e.total_uti_patients  = numel(patTotUTI);
    e.valid_uti_patients  = numel(patValUTI);

    e.total_uti_culture_ids    = cultTot;
    e.valid_uti_culture_ids    = cultVal;
    e.discarded_uti_culture_ids = setdiff(cultTot,cultVal);

    e.total_uti_cultures   = numel(cultTot);
    e.valid_uti_cultures   = numel(cultVal);
    e.discarded_uti_cultures = e.total_uti_cultures - e.valid_uti_cultures;

    S.exp.(thisExp) = e;

    % running isolate totals
    S.total_fobt_isolates = S.total_fobt_isolates + e.total_fit_isolates;
    S.total_uti_isolates = S.total_uti_isolates + e.total_uti_isolates;
    S.valid_fobt_isolates = S.valid_fobt_isolates + e.valid_fit_isolates;
    S.valid_uti_isolates = S.valid_uti_isolates + e.valid_uti_isolates;

    % build grand-patient sets
    fobtPatientSet = [fobtPatientSet ; patTotFIT];
    utiPatientSet = [utiPatientSet ; patTotUTI];    
    fobtPatientSetV = [fobtPatientSetV ; patValFIT];
    utiPatientSetV = [utiPatientSetV ; patValUTI];
    utiCultureSet = [utiCultureSet ; cultTot];
    utiCultureSetV = [utiCultureSetV ; cultVal];
end


% ------------------------------------------------------------
% 6. GRAND-TOTAL PATIENT COUNTS (across all experiments)

S.n_fobt_patients_total = numel(unique(fobtPatientSet));
S.fobt_patients_total = unique(fobtPatientSet);

S.n_uti_patients_total = numel(unique(utiPatientSet));
S.uti_patients_total = unique(utiPatientSet);

S.n_patients_total     = numel(unique([fobtPatientSet ; utiPatientSet]));
S.patients_total     = unique([fobtPatientSet ; utiPatientSet]);

S.n_uti_cultures_total = numel(unique(utiCultureSet));
S.uti_cultures_total = unique(utiCultureSet);

S.n_fobt_patients_valid = numel(unique(fobtPatientSetV));
S.fobt_patients_valid = unique(fobtPatientSetV);

S.n_uti_patients_valid = numel(unique(utiPatientSetV));
S.uti_patients_valid = unique(utiPatientSetV);

S.n_uti_cultures_valid = numel(unique(utiCultureSetV));
S.uti_cultures_valid = unique(utiCultureSetV);
end