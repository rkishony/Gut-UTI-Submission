%% Reset

clear all
close all
clc
add_required_paths
folders = get_folders();
record_value([]);


%% Code flow configuration
global CONFIG  %#ok<GVMIS>

% Control whether to recreate data from raw, or load from mat:
force_create = true;

% Control figure creation and printing:
CONFIG = struct;
CONFIG.create_figures = true;
CONFIG.print_figures = true;  % only active when create_figures=true

% Control creation of xls file with all numeric values:
CONFIG.save_values_to_xls = true;


%% Load all data

load_or_create('load_resfinder', 'resfinder.mat', force_create)
load_or_create('load_allMaccabiData', 'allMaccabiData.mat', force_create)
samepatcovbreadth = load_alignment(true, false);

[PHEN1, PHEN1_Abs] = load_CommunityResPhen('CommunityPhenotyping');
[PHEN2, PHEN2_Abs] = load_CommunityResPhen('CommunityPhenotypingB');

%% Set Parameters

% Time spans for linking events
Timespan = struct;
Timespan.BAC_FOBT_inclusion     = [ -1, 365];
Timespan.DIAG_BAC_coinciding    = [-14,  14];
Timespan.BAC_FOBT_coinciding    = [ -3,   3];

% MHS resistance measurements
BAC_or_SMP = 'SMP';  % Associate FOBT with bac res (BAC), or with uti sample (SMP)
I_as_RS    = 'R';    % Whether "I" resistance should be translated to "S" or "R"

% Meta-genomics
Etha0_th                        = 1e-4;  % Meta-Genomic FOBT-UTI genome coverage threshold 
                                         % (to define 'identity')
% Resistance thresholds
Res_th = struct;
Res_th.cov                      = 0.05;  % Genome-normalized gene coverage
Res_th.mut                      = 1.5;   % Number of resistant mutations (the fraction doesn't matter, just for visualization)
Res_th.phen                     = 0.05;  % Fraction of resistant colonies

% Meta-Genomic Resfinder
Genome_cov_limit                = 1;     % Limit for discarding samples with too-low coverage

RF_params = struct;
RF_params.presumed_genome_len   = 5e6;   % To estimate genome coverage
RF_params.coverage_thr          = 60;    % Discard Resfinder hits with lower coverage
RF_params.identity_thr          = 95;    % Discard Resfinder hits with lower identity score
RF_params.exclude_singleton     = false; % Exclude Resfinder files with single row (legacy)

% Meta-Phenotyping
PHEN_params = struct;
PHEN_params.max_colony_count    = 60;
PHEN_params.min_colony_count    = 25;

% Colors of sen/res:
CLRS.sen = [0.6144 0.5874 0.7379];
CLRS.res = [0.8415 0.5887 0.5107];

fmt = struct;
fmt.RF_params.presumed_genome_len = "%1.0e";
fmt.Etha0_th = "%1.0e";
record_value('Parameters', 'params', var2struct(Res_th, RF_params, PHEN_params, Etha0_th), fmt)

%% Process PHEN
% Normalize the meta-resistance measurements by the MC control,
% and remove faulty measurements

[PHEN1, PHEN1_Entities] = process_PHEN(PHEN1, PHEN1_Abs, ...
    PHEN_params.max_colony_count, PHEN_params.min_colony_count);
[PHEN2, PHEN2_Entities] = process_PHEN(PHEN2, PHEN2_Abs, ...
    PHEN_params.max_colony_count, PHEN_params.min_colony_count);

[PHEN, PHEN_Entities] = combine_phens(PHEN1, PHEN2, PHEN1_Entities, PHEN2_Entities);

%% Define our patient cohort (POP) and add cross-table variables

merge_FOBT_with_MHS_memorize = file_memorize(@merge_FOBT_with_MHS, force_create);

[POP, FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN] = ...
    merge_FOBT_with_MHS_memorize(FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN, Timespan);

if strcmp(BAC_or_SMP, 'BAC')
    PUC = BAC;
else
    PUC = SMP;
end


%% Process RESFDR
% We pivot Resfinder's lists of coverage (RESFDR_COV) and mutations (RESFDR_MUT)
% into matrices of MAX coverages and SUM mutations, 
% for Alleles, Genes, and Drug (which we collectively term 'entity') for 
% each sample.

% Keep only resfinder data matching fobt:
RF = RESFDR(RESFDR.jFOBT>0,:);

% Add label and genome coverage for each sample:
RF.Label = arrayfun(@(x) sprintf('F%d',x), RF.FOBT_ID, 'UniformOutput',false);
RF.avgcov = RF.bpTotal / RF_params.presumed_genome_len;
RF.isGood = RF.avgcov > Genome_cov_limit;  % Define good (valid) samples

% MAX-pivot RESFDR_COV to create cov_Allele, cov_Gene, cov_Drug (new cols in RF):
ok = ~(RESFDR_COV.Singleton & RF_params.exclude_singleton);
ok = ok & RESFDR_COV.Template_Coverage > RF_params.coverage_thr;
ok = ok & RESFDR_COV.Query_Identity > RF_params.identity_thr;
[RF, RF_Entities_cov] = pivot_resfinder(RF, RESFDR_COV(ok,:), @max, 'cov_', RF.avgcov);

% SUM-pivot RESFDR_MUT to create mut_Allele, mut_Gene, mut_Drug (new cols in RF):
ok = ~(RESFDR_MUT.Singleton & RF_params.exclude_singleton);
[RF, RF_Entities_mut] = pivot_resfinder(RF, RESFDR_MUT(ok,:), @sum, 'mut_');

% RF_Entities is a struct with field names identical to the data columns of
% resfinder's RF, listing for each column the specific entities (drugs,
% genes, or alleles) that are provided in the RF column.
RF_Entities = joinStruct(RF_Entities_cov, RF_Entities_mut);

% Define TMP+SMX combination (AND gate, namely 'min' coverage):
j_tmp_smx = ismember(RF_Entities.cov_Drug.Drug, {'trimethoprim','sulphonamide'});
RF.cov_TMP_SMX = min(RF.cov_Drug(:,j_tmp_smx),[],2);

clear ok RF_Entities_cov RF_Entities_mut j_tmp_smx


%% Define Cohorts
% We can now list our patient cohort with their demographics and 
% which lab data we have for each.

% Cohorts is a table same height as POP.

Cohorts = table;

% Base
Cohorts.All = POP.anyPUC & POP.FOBTisFirst;
assert(all(Cohorts.All))

% Coinciding
Cohorts.CollectedAsCoinciding = POP.FOBT_ID<1000;
Cohorts.CollectedAsNonCoinciding = POP.FOBT_ID>=1000;
Cohorts.Coinciding = POP.CoincidingPUC;
Cohorts.NonCoinciding = ~POP.CoincidingPUC;

% Demographic risk group
Cohorts.F5080 = inrange(POP.age, [50 80]) & POP.Gender=='F';

% Patients with Resfinder data
Cohorts.isRF = POP.jRF > 0;
Cohorts.isGoodRF = safe_index(RF.isGood, POP.jRF) > 0;
Cohorts.isGoodRF_F5080 = Cohorts.isGoodRF & Cohorts.F5080;
Cohorts.isGoodRF_NonCo = Cohorts.isGoodRF & Cohorts.NonCoinciding;
Cohorts.isGoodRF_CollectedAsNonCo = Cohorts.isGoodRF & Cohorts.CollectedAsNonCoinciding;

% Patients with Metaphenotyping data
Cohorts.isPHEN = POP.jPHEN > 0;

% Patients with UTI isolates
sumUTIiso = countMatches(POP, UTIpop, 'FOBT_ID');
Cohorts.isUTIiso = sumUTIiso > 0;

% Check RF-POP consitency:
assert(seteq(find(RF.isGood & RF.FOBTisFirst & ~RF.isrep), POP.jRF(Cohorts.isGoodRF)))


%% Output/record cohort sizes:

% Cohort / POP
countCohorts = structfun(@sum, table2struct(Cohorts, 'ToScalar', true), 'UniformOutput',false);
countCohorts.women = sum(POP.Gender(Cohorts.All)=='F');
countCohorts.men = sum(POP.Gender(Cohorts.All)=='M');
countCohorts.above_50 = sum(POP.age>=50);
fprintf('\n')
record_value('Cohort count (POP)', [], countCohorts)

% FOBT
countFOBT = struct;
ok = FOBT.FOBTisFirst;
countFOBT.AllPatients = sum(ok);

countFOBT.NoFollowingPUC = sum(ok & ~FOBT.anyPUC);

coinc = ok & FOBT.FOBT_ID<1000;
countFOBT.CollectedAsCoinciding = sum(coinc);
countFOBT.CoincidingStartDate = datestr(min(FOBT.FOBT_Date(coinc)), 'yyyy-mm-dd');
countFOBT.CoincidingEndDate = datestr(max(FOBT.FOBT_Date(coinc)), 'yyyy-mm-dd');
countFOBT.CollectedAsCoincidingAndIsCoinciding = sum(coinc & FOBT.CoincidingPUC);

nonco = ok  & FOBT.FOBT_ID>=1000;
countFOBT.CollectedAsNonCoinciding = sum(nonco);
t = FOBT.FOBT_Date(nonco);
t0 = datenum('2021-01-01');  % any date between the collection periods.
countFOBT.NonCoincidingStartDate1 = datestr(min(t), 'yyyy-mm-dd');
countFOBT.NonCoincidingEndDate1 = datestr(max(t(t<t0)), 'yyyy-mm-dd');
countFOBT.NonCoincidingStartDate2 = datestr(min(t(t>=t0)), 'yyyy-mm-dd');
countFOBT.NonCoincidingEndDate2 = datestr(max(t), 'yyyy-mm-dd');

countFOBT.CollectedAsNonCoincidingButCoinciding = sum(nonco & FOBT.CoincidingPUC);

record_value('Feacal samples (FOBT)', [], countFOBT)

% SMP
countSMP = struct;
ok = SMP.FOBTisFirst;
countSMP.LastRecordedPUCDate = datestr(max(SMP.SampleDate(ok)), 'yyyy-mm-dd');
ok = ok & inrange(SMP.DateDiff, Timespan.BAC_FOBT_inclusion);
countSMP.TotalSMPwithinOneYear = sum(ok);
countSMP.TotalBACwithinOneYear = sum(SMP.nBac(ok));

record_value('MHS positive urine cultures (SMP)', [], countSMP)

% PHEN_lables (labeling each line of PHEN)
PHEN_labels = table;
PHEN_labels.all = true(height(PHEN), 1);
PHEN_labels.ismergedMC = PHEN.ismerged(:,1);
PHEN_labels.below_threshold_MC = PHEN.clean_counts(:,1) < PHEN_params.min_colony_count;
PHEN_labels.isrep = PHEN.isrep;
PHEN_labels.inPOP = ismember(PHEN.FOBT_ID, POP.FOBT_ID);
PHEN_labels.inPOP_notrep = PHEN_labels.inPOP & ~PHEN_labels.isrep;
PHEN_labels.inPOP_notrep_notismergedMC = PHEN_labels.inPOP_notrep & ~PHEN_labels.ismergedMC;
PHEN_labels.inPOP_notrep_notismergedMC_lowcount = PHEN_labels.inPOP_notrep_notismergedMC & PHEN_labels.below_threshold_MC;
PHEN_labels.inPOP_ok = PHEN_labels.inPOP_notrep_notismergedMC & ~PHEN_labels.below_threshold_MC;

countPHEN = structfun(@sum, table2struct(PHEN_labels, 'ToScalar', true), 'UniformOutput',false);

record_value('Community phenotyping (PHEN)', [], countPHEN)

% UTIpop
UTIpop_labels = table;
UTIpop_labels.all = UTIpop.FOBTisFirst | isnan(UTIpop.FOBT_ID);
assert(all(UTIpop_labels.all))  % Because we did not collect from the patients with multiple FOBTs
UTIpop_labels.WithFOBT = UTIpop.FOBTisFirst;
UTIpop_labels.InPOP = UTIpop.FOBTisFirst & ismember(UTIpop.FOBT_ID, POP.FOBT_ID);
UTIpop_labels.InPOPandWithinYear = UTIpop_labels.InPOP & inrange(UTIpop.DateDiff, Timespan.BAC_FOBT_inclusion);

countUTI = structfun(@sum, table2struct(UTIpop_labels, 'ToScalar', true), 'UniformOutput',false);

record_value('Collected UTI cultures (UTIpop)', [], countUTI)


%% Construct matching FOBT-PUC data
%
% Unified Genotypic (RF) and phenotypic (PHEN) analysis of FOBT - PUC
% matching.
%
% fobt_val  [nF, nD]           FOBT values (RF/PHEN) by faecal sample and drug
% v123count [nF, nD, nT, 3]    PUCs are organized as array of MHS counts 
%                              per faecal sample, drug, and time bin post
%                              sample, for 
%                              1: Null (PUC without measurement)
%                              2: Sen (Sensitive PUC)
%                              3: Res (Resistant PUC)
% 
% This organization is done per Population Condition, which specifies which
% FOBT and PUC data to include.


% ##### Define MHS resistance groups #####

% Define sets of HMS-measured drugs (sets of RES_All's indices), which
% can then be used to define PUC resistance to be matched with genotypic or 
% phenotypic FOBT resistances.

MHS_ResGroups = struct;

% Define HMS-measured drug groups for each Resfinder resistance
drugs = unique(RES_All.ResFinder);
for i = 1:numel(drugs)
    MHS_ResGroups.(drugs{i}) = find(strcmp(drugs{i}, RES_All.ResFinder));
end

% Add commonly measured drugs as individual groups
for i = 1:20 
    MHS_ResGroups.(makeValidName(RES_All.Name{i})) = i;
end

% Add 'CEF' group:
[~, MHS_ResGroups.CEF] = ismember({'Cefazolin', 'Ceftazidime', 'Cefuroxime axetil'}, RES_All.Name);
% 'Cephalexin'  maybe used in community

clear drugs i


% ##### Define populations #####

puc_base = @(puc) true(height(puc),1); % all
puc_diag = @(puc) puc.isuti;
puc_diff_sp = @(puc) ~puc.isSameSpeciesAsCoinciding;

PopulationConditions = cell2table({
  %
  %                                              PUC              time
  % tag     type   POP mask         PUC mask     time method*     frame
  % -----------------------------------------------------------------------
    "base"  "RF"   "isGoodRF"       puc_base     "all"            [0 365]
    "F5080" "RF"   "isGoodRF_F5080" puc_base     "all"            [0 365]
    "NonCo" "RF"   "isGoodRF_NonCo" puc_base     "all"            [0 365]
    "ClctNonCo" "RF"   "isGoodRF_CollectedAsNonCo" ...
                                    puc_base     "all"            [0 365]
    "diagU" "RF"   "isGoodRF"       puc_diag     "all"            [0 365]
    "difSp" "RF"   "isGoodRF"       puc_diff_sp  "all"            [0 365]
    "first" "RF"   "isGoodRF"       puc_base     "first_measured" [0 365]
    "base"  "PHEN" "isPHEN"         puc_base     "all"            [0 365]
    "base7"  "PHEN" "isPHEN"         puc_base     "all"            [7 365]
    "base21"  "PHEN" "isPHEN"         puc_base     "all"            [21 365]
    }, ...
    'VariableNames', {'tag', 'type', 'pop_mask', 'puc_mask', 'puc_time_method', 'timeframe'});
nCond = height(PopulationConditions);

%  * Options listed in: `help count_UTI_resistances`
%    This field is indicated for subsequent analysis (not used here)


% ##### Define mapping from FOBT's ResFinder/PHEN names to MHS measurements #####

DrugMapping = cell2table({ ...
   %
   %                        =========== FOBT Data ===========  === PUC data ===
   %                        RF/PHEN          RF/PHEN           Indicate field in    Resistance
   % type   Pretty          Entity           Field             MHS_ResGroups        type
   % --------------------------------------------------------------------------------------------------------
    "RF"   "\beta-lactam"  "cov_Drug"       "beta_lactam"     "beta_lactam"         "cov"
    "RF"   "FOS"           "cov_Drug"       "fosfomycin"      "fosfomycin"          "cov"
    "RF"   "TMP+SMX"       "cov_TMP_SMX"    ""                "Trimethoprim_Sulfa"  "cov"
    "RF"   "CPR (gene)"    "cov_Drug"       "quinolone"       "quinolone"           "cov"
    "RF"   "CPR (mut)"     "mut_Drug"       ""                "quinolone"           "mut"
    "RF-"  "AGs"           "cov_Drug"       "aminoglycoside"  "aminoglycoside"      "cov"  % Skip. Not enough data

    "PHEN" "CEF"           "clean_counts_n" "CEF_A"           "CEF"                 "phen"
    "PHEN" "NIT"           "clean_counts_n" "NIT-8_B"         "Nitrofurantoin"      "phen"
    "PHEN" "FOS"           "clean_counts_n" "FOS-16_B"        "fosfomycin"          "phen"
    "PHEN" "TMP"           "clean_counts_n" "TMP_A"           "Trimethoprim_Sulfa"  "phen"
    "PHEN" "CPR"           "clean_counts_n" "CPR_A"           "quinolone"           "phen"
    }, ...
    'VariableNames', {'type', 'pretty', 'entity', 'field', 'ab_group', 'res_type'});


% Define time dins for PUC-FOBT date diff:
DayEdges = [-inf, -1.5:365.5, inf];

% Initiate per-condicion collected data:
CondResults = struct;

% ##### Calculate patient-paired PUC x FindRes/PHEN #####

for iCond = 1:nCond
    cond_row = PopulationConditions(iCond,:);
    cond_results = build_matching_fobt_puc_data(DrugMapping, cond_row, ...
        POP, Cohorts, PUC, MHS_ResGroups, I_as_RS, DayEdges, RF, RF_Entities, PHEN, PHEN_Entities);
    CondResults.(cond_results.label) = cond_results;
end

clear iCond cond_row cond_results

%%
% ==============================================================================================
% ======================================= Start Figures ========================================
% ==============================================================================================

%% Plot CDF

% Set FOBT thresholds and plot parameters

% Example FOBT to mark on CDF of PHEN.TMP
% (Must match the images created by plot_examples.m on the Dropbox folder)

% Define FOBT IDs to be marked in PHEN_TMP to indicate colony image examples
FOBT_IDS_Example = [
  % Sen     Res
    240     407
    467     183
    133     262
    ];
markers = '^so';

res_th_type2params = cell2table({...
  %  res_type   res_th       x_lim       is_log
    'cov'       Res_th.cov   [1e-4 1]    1
    'mut'       Res_th.mut   [0    4]    0
    'phen'      Res_th.phen  [1e-2 1]    1
},'VariableNames',{'res_type', 'res_th', 'x_lim', 'islog'}, 'RowNames',{'cov', 'mut', 'phen'});

% Plot FOBT CDF by PUC resistance
cond_labels = string(fields(CondResults));
for iCond = 1:numel(cond_labels)
    cond_label = cond_labels(iCond);
    cond_results = CondResults.(cond_label);
    ConfMatDrugs = nan(2,2,cond_results.nDrugs);

    fprintf('\n%25s %s\n', '"' + cond_label + '"', repmat('-',1,58));
    nDrugs = cond_results.nDrugs;
    fobt_ths = nan(1, nDrugs);
    for ijDrugs = 1:nDrugs
        jDrug = cond_results.jDrugs(ijDrugs);
        dmap = table2struct(DrugMapping(jDrug,:));
        param = res_th_type2params(dmap.res_type,:);
        label_drug = strjoin([cond_label, dmap.pretty], '_');
        if dmap.pretty == "\beta-lactam"
            nres = 3;
        else
            nres = 1;
        end
        fobt_th = param.res_th;
        fobt_ths(ijDrugs) = fobt_th;

        fprintf('%25s', dmap.pretty)

        [num_puc, num_puc_measured, num_puc_res] = count_UTI_resistances(...
            cond_results.v123count(:, ijDrugs, :, :), DayEdges, cond_results.timeframe, cond_results.puc_time_method);
        fobt_val = cond_results.fobt_val(:,ijDrugs);

        % Get the confusion matrix and p-values for association betwebe the
        % FOBT value and the MHS recorded resistance:
        %%% TODO: perhaps good to unify `calc_confmat_over_drugs_and_time_with_permutations`
        % and `get_fobt_puc_confusion_mat_and_association_statistics`

        [conf_mats, chosen, Fisher_p_value, KS_p_value] = ...
            get_fobt_puc_confusion_mat_and_association_statistics(...
            fobt_val, fobt_th, num_puc_measured, num_puc_res, nres);
        fprintf(' Fisher=%-24s KS=%-24s ', sprintf('%7.5f ', Fisher_p_value), sprintf('%7.5f ', KS_p_value))
        ConfMatDrugs(:, :, ijDrugs) = conf_mats(:,:,1);

        % Plot the cdf of FOBT values for patients with sen / res PUCs:
        if any(sum(chosen,1)<4)
            fprintf('  (Too little counts, skip cdf)')
        elseif any(cond_label == ["RF_base" "PHEN_base"])
            fig_num = 10000 + 3000 + iCond*100 + ijDrugs;
            [cf, h, X, Y] = fig_cdf_res_sen(fig_num, cond_results.type, fobt_val, chosen, dmap.pretty, fobt_th, param.islog, param.x_lim, CLRS.sen, CLRS.res);
            if ~isempty(h)
                if ~ismember(label_drug, ["RF_base_\beta-lactam", "RF_base_TMP+SMX", "PHEN_base_TMP", "PHEN_base_CPR"])
                    delete(h.legend)
                end
                if ismember(label_drug, ["RF_base_\beta-lactam", "PHEN_base_TMP", "PHEN_base_CPR"])
                    h.val_th_text.Position(2) = 0.94;
                end
                if label_drug == "PHEN_base_TMP"
                    for j_sr = 1:2
                        if j_sr == 1, clr = CLRS.sen; else clr = CLRS.res; end
                        for ke = 1:size(FOBT_IDS_Example,1)
                            j_d = find(PHEN.FOBT_ID(cond_results.jFDATA) == FOBT_IDS_Example(ke,j_sr));
                            plot(X(j_d, j_sr), Y(j_d, j_sr), markers(ke), 'MarkerEdgeColor', clr * 0.7, 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.5)
                        end
                    end
                end
                if 0
                    for j_d = 1:size(X,1)
                        id = PHEN.FOBT_ID(cond_results.jFDATA(j_d));
                        id = [' ' num2str(id)];
                        if ~isnan(X(j_d,1))
                            text(X(j_d,1), Y(j_d,1), id, 'color','b')
                        elseif ~isnan(X(j_d,2))
                            text(X(j_d,2), Y(j_d,2), id, 'color','r')
                        end
                    end
                end
                print_figure(cf, sprintf('fig_cdf_res_sen_%s_%s', cond_label, strrep(dmap.pretty, '\', '')));
            end
        end
        fprintf('\n')
    end
    CondResults.(cond_label).fobt_ths = fobt_ths;
    CondResults.(cond_label).ConfMatDrugs = ConfMatDrugs;
end
fprintf('\n')

clear ConfMatDrugs



%% Barcharts of PUC resistance risk stratified by genotypic/phenotypic FOBT resistance
% Based on the confusion matrices for each condition and each drug (CondResults.(cnd_tag).ConfMatDrugs), 
% plot a barchat of risk of PUC resistance for sensitive/resistant FOBT.

% Choose rows from PopulationConditions table to plot as bars.
% (combine two rows to create a barchart with two types of measurements (like for t1>0, t1>30)
CondLabelGroups = {"RF_base", "PHEN_base", "PHEN_base7", "PHEN_base21"};
y_lim  = [1.1, 1, 1, 1];

for i = 1:numel(CondLabelGroups)
    cond_lbls = CondLabelGroups{i};
    cond_arr = array_struct_fields(CondResults, cond_lbls);
    label = strjoin([cond_arr.label],'_');
    rf_or_phen = cond_arr(1).type;
    assert(all(strcmp(rf_or_phen, [cond_arr.type])))
    drugs = DrugMapping.pretty(DrugMapping.type == rf_or_phen);
    drugs = cellstr(drugs);
    for k = 1:numel(drugs)
        f = find(drugs{k} == '(', 1);
        if ~isempty(f)
            drugs{k} = {drugs{k}(1:f-1), drugs{k}(f:end)};
        end
    end

    if numel(cond_lbls) == 2
        assert(isequal(cond_arr(1).Cond(:,[2 4:end]), cond_arr(2).Cond(:,[2 4:end])))
        t1_vec = [cond_arr(1).Cond.timeframe(1), cond_arr(2).Cond.timeframe(1)];
        t1_vec(t1_vec==-1) = 0;  % account for 1-off clinical tolerance
        t1_labels = arrayfun(@(t) sprintf('\\Deltat\\geq%d', t), t1_vec, 'UniformOutput',false);
        labels = {{}, t1_labels, drugs};
    else
        labels = {{}, {}, drugs};
    end
    conf_mat = cat(4, cond_arr.ConfMatDrugs);
    conf_mat = permute(conf_mat,[1, 2, 4, 3]);

    [cf, h] = fig_bar_res_sen(10200+i, conf_mat, labels, rf_or_phen, CLRS.sen, CLRS.res, y_lim(i));
    if ~isempty(h)
        if ~any(label == ["RF_first" "RF_base" "PHEN_base"])
            % no need. will have same legend from other panels.
            delete(h.legend.text)  
            delete(h.legend.graphics)  
        end
    end
    print_figure(cf, sprintf('fig_bar_res_sen_%s', label));
end


%% FOBT-UTI drug-combined odds ratio over time
rng(54)

% Choose conditions from PopulationConditions:
cond_and_legend = {
    "RF_base"   "All"                       1   1
    "RF_F5080"  "Female 50-80 year old"     0   1
    "RF_NonCo"  "Non-coinciding cohort"     0   1
    "RF_diagU"  "UTI Diagnostic"            0   1
    "RF_difSp"  "Different species"         0   1
    "RF_first"  "First PUC"                 0   1
    "PHEN_base" "All-Pheno"                 1   0
    };
randomization = 'complete';  % 'complete' / 'restricted'
OR_v = [];
OR_bs_prctile_v = [];
OR_perm_P_val = [];
OR_perm_prctile_v = [];
for iCond = 1:size(cond_and_legend,1)
    cond_lbl = cond_and_legend{iCond,1};
    cond_results = CondResults.(cond_lbl);
    
    % Choose which drugs to use for aggregated odds ratio:
    if cond_results.type == "RF"
        % Remove CPR(gene):
        jd = find(DrugMapping.pretty(cond_results.jDrugs) ~= "CPR (gene)");
    else
        % All drugs:
        jd = 1:cond_results.nDrugs;
    end

    % FOBT Data   [nF, nD]
    F3 = convert_fobt_val_to_F3(cond_results.fobt_val, cond_results.fobt_ths);
    F3 = F3(:,jd);

    % Get PUC counts
    v123count = cond_results.v123count;  % [nF, nD, nT, 3]  1 NULL, 2 SEN, 3 RES.
    v123count = v123count(:,jd,:,:);

    % Choose errobar percentile
    prctiles = [16 50 84];  % One std
    
    % Remove the flanking inf timebins (-inf, +inf):
    assert(DayEdges(1)==-inf & DayEdges(end)==inf)
    day_edges = DayEdges(2:end-1);
    v123count = v123count(:,:,2:end-1,:);  % remove inf bins
    days = (day_edges(1:end-1) + day_edges(2:end)) / 2;

    % Get V4: (None, Null, Sen, Res for each sample x drug x time)
    V4 = convert_v123count_to_v4(v123count, [], cond_results.puc_time_method, 3);

    % Choose days:
    jT = find(days==cond_results.timeframe(1)):56:365;
    days = days(jT);
    V4 = V4(:,:,jT);

    % Calculate nominal and boostrap and permuted odds ratios (OR):
    [CM, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile, fit_nominal, fit_bs_prctile] = ...
        calc_F3_V4_association_across_time_and_drugs(F3, V4, days, prctiles, randomization);
    
    % Verify consistency with the t=0 with the calc of association without time:
    assert(isequal(CM(2:3,3:4,:,1), CondResults.(cond_lbl).ConfMatDrugs(:,:,jd)))

    % Store the first time point
    OR_v(:,:,iCond) = OR;
    OR_bs_prctile_v(:,:,iCond) = OR_bs_prctile;
    OR_perm_prctile_v(:,:,iCond) = OR_perm_prctile;
    OR_perm_P_val_v(:,:,iCond) = OR_perm_P_val;
    if cond_and_legend{iCond,3}
        cf = fig_oddsratio_vs_time(30+iCond, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile, days);
        print_figure(cf, char(strjoin(['fig_oddsratio_vs_time' cond_lbl],'_')));
    end
end


%% Drug-combined Odds Ratio
ok = cell2mat(cond_and_legend(:,4)) == 1;
cf = fig_oddsratio_vs_time_all(30, OR_v(:,:,ok), OR_perm_P_val_v(:,:,ok), OR_perm_prctile_v(:,:,ok), OR_bs_prctile_v(:,:,ok), days, string(cond_and_legend(ok,2)));
print_figure(cf, 'fig_fig_oddsratio_vs_time_all')


%% Heatmap of resistance gene coverage by sample

ok = RF.isGood & RF.FOBTisFirst & ~RF.isrep;
[cf, gene_cov_map.num_samples_not_plotted] = fig_res_gene_coverage_matrix(101, RF(ok,:), RF_Entities, 'Gene', 'Drug', Res_th.cov);
record_value('Gene coverage heatmap', [], gene_cov_map)
print_figure(cf, 'fig_res_gene_coverage_matrix', [], [], 'Arial')


%% Patient statistics histograms

ValidSMP = SMP(inrange(SMP.DateDiff, Timespan.BAC_FOBT_inclusion) & ...
    ismember(SMP.FOBT_ID,POP.FOBT_ID), :);
cf = fig_pop_stat(102, POP, ValidSMP, Cohorts.Coinciding, Cohorts.NonCoinciding);
print_figure(cf, 'fig_pop_stat');


%% PUCs timelines per patient, indicating genomic matching

cf = fig_PUC_genomic_match_timelines(103, UTIpop, SMP, samepatcovbreadth, Etha0_th);
print_figure(cf, 'fig_PUC_genomic_match_timelines')


%% Resfinder consitency in same-sample replicates versus different samples

cRF = RF(RF.isGood & RF.FOBTisFirst, :);
[rf_reps.tot_num_reps, rf_reps.tot_unique_ids, dists, rf_reps.p_vals] = ...
    analyze_rep_vs_non_rep(cRF.cov_Drug > Res_th.cov, cRF.FOBT_ID);
rf_reps.avgs = cellfun(@mean, dists);
record_value('Resfinder replicates', [], rf_reps);

cf = fig_FIT_gen_res_profile_replicates(104, dists, rf_reps.p_vals);

print_figure(cf, 'fig_FIT_gen_res_profile_replicates')


%% Create Excel with recorded values

values = record_value();
if CONFIG.save_values_to_xls
    writetable(values, 'values.xls')
end