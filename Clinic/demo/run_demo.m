function run_demo()
% RUN_DEMO  End-to-end demonstration of the Clinic analysis on synthetic data.
%
% This runs the *real* Clinic analysis functions (merge_FOBT_with_MHS,
% build_matching_fobt_puc_data, the confusion-matrix / odds-ratio statistics
% and the figure functions) on a small SYNTHETIC dataset produced by
% make_synthetic_data.m. No protected patient data are used.
%
% Outputs are written to Clinic/demo/output/:
%   figures/*.svg   - the generated figures
%   values.xls      - the recorded numeric values
%
% Usage (from MATLAB):
%   cd Clinic/demo
%   run_demo
%
% Expected run time: about 1-2 minutes on a normal desktop (most of it in the
% permutation/bootstrap odds-ratio analysis).

demo_dir = fileparts(mfilename('fullpath'));
clinic_root = fileparts(demo_dir);

% ---- Paths --------------------------------------------------------------
addpath(genpath(fullfile(clinic_root, 'utils')));
addpath(fullfile(clinic_root, 'analysis_funcs'));
addpath(fullfile(clinic_root, 'loading_scripts'));
addpath(fullfile(clinic_root, 'fig_scripts'));
addpath(demo_dir);

% ---- Output location (print_figure writes to ./figures) ----------------
out_dir = fullfile(demo_dir, 'output');
if ~exist(fullfile(out_dir, 'figures'), 'dir')
    mkdir(fullfile(out_dir, 'figures'));
end
start_dir = pwd;
cleanupObj = onCleanup(@() cd(start_dir));
cd(out_dir);

% ---- Config -------------------------------------------------------------
global CONFIG %#ok<GVMIS>
CONFIG = struct('create_figures', true, 'print_figures', true, 'save_values_to_xls', true);
record_value([]);

% ---- Load synthetic data (create it if missing) ------------------------
data_file = fullfile(demo_dir, 'demo_data.mat');
if ~exist(data_file, 'file')
    make_synthetic_data(data_file);
end
L = load(data_file);
FOBT=L.FOBT; BAC=L.BAC; SMP=L.SMP; IDS=L.IDS; DIAG=L.DIAG; UTIpop=L.UTIpop; %#ok<NASGU>
RES_All=L.RES_All; RESFDR=L.RESFDR; RESFDR_COV=L.RESFDR_COV; RESFDR_MUT=L.RESFDR_MUT;
PHEN=L.PHEN; PHEN_Entities=L.PHEN_Entities;

% ---- Parameters (same values as generateFigures.m) ---------------------
Timespan = struct;
Timespan.BAC_FOBT_inclusion  = [-1, 365];
Timespan.DIAG_BAC_coinciding = [-14, 14];
Timespan.BAC_FOBT_coinciding = [-3, 3];

I_as_RS = 'R';

Res_th = struct('cov',0.05, 'mut',1.5, 'phen',0.05);
RF_params = struct('presumed_genome_len',5e6, 'coverage_thr',60, ...
    'identity_thr',95, 'exclude_singleton',false);
Genome_cov_limit = 1;
PHEN_params = struct('max_colony_count',60, 'min_colony_count',25);

CLRS.sen = [0.6144 0.5874 0.7379];
CLRS.res = [0.8415 0.5887 0.5107];

record_value('Parameters', 'params', var2struct(Res_th, RF_params, PHEN_params));

% ---- Process ResFinder (pivot coverage & mutations) --------------------
RF = RESFDR;
RF.Label = arrayfun(@(x) sprintf('F%d',x), RF.FOBT_ID, 'UniformOutput',false);
RF.avgcov = RF.bpTotal / RF_params.presumed_genome_len;
RF.isGood = RF.avgcov > Genome_cov_limit;

ok = ~(RESFDR_COV.Singleton & RF_params.exclude_singleton);
ok = ok & RESFDR_COV.Template_Coverage > RF_params.coverage_thr;
ok = ok & RESFDR_COV.Query_Identity > RF_params.identity_thr;
[RF, RF_Entities_cov] = pivot_resfinder(RF, RESFDR_COV(ok,:), @max, 'cov_', RF.avgcov);

ok = ~(RESFDR_MUT.Singleton & RF_params.exclude_singleton);
[RF, RF_Entities_mut] = pivot_resfinder(RF, RESFDR_MUT(ok,:), @sum, 'mut_');
RF_Entities = joinStruct(RF_Entities_cov, RF_Entities_mut);

j_tmp_smx = ismember(RF_Entities.cov_Drug.Drug, {'trimethoprim','sulphonamide'});
RF.cov_TMP_SMX = min(RF.cov_Drug(:,j_tmp_smx),[],2);
RESFDR = RF;   % keep the enriched table
clear ok RF_Entities_cov RF_Entities_mut j_tmp_smx

% ---- Merge FOBT with MHS-style tables ----------------------------------
[POP, FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN] = ...
    merge_FOBT_with_MHS(FOBT, BAC, SMP, IDS, DIAG, RESFDR, PHEN, Timespan); %#ok<ASGLU>
RF = RESFDR(RESFDR.jFOBT>0,:);
PUC = SMP;

% ---- Cohorts ------------------------------------------------------------
Cohorts = table;
Cohorts.All           = POP.anyPUC & POP.FOBTisFirst;
Cohorts.Coinciding    = POP.CoincidingPUC;
Cohorts.NonCoinciding = ~POP.CoincidingPUC;
Cohorts.F5080         = inrange(POP.age,[50 80]) & POP.Gender=='F';
Cohorts.isGoodRF      = safe_index(RF.isGood, POP.jRF) > 0;
Cohorts.isGoodRF_F5080= Cohorts.isGoodRF & Cohorts.F5080;
Cohorts.isPHEN        = POP.jPHEN > 0;

counts = structfun(@sum, table2struct(Cohorts,'ToScalar',true), 'UniformOutput',false);
counts.women = sum(POP.Gender(Cohorts.All)=='F');
counts.men   = sum(POP.Gender(Cohorts.All)=='M');
record_value('Cohort count (POP)', [], counts)

% ---- MHS resistance groups ---------------------------------------------
MHS_ResGroups = struct;
drugs = unique(RES_All.ResFinder);
for i = 1:numel(drugs)
    MHS_ResGroups.(makeValidName(drugs{i})) = find(strcmp(drugs{i}, RES_All.ResFinder));
end
for i = 1:20
    MHS_ResGroups.(makeValidName(RES_All.Name{i})) = i;
end
[~, MHS_ResGroups.CEF] = ismember({'Cefazolin','Ceftazidime','Cefuroxime axetil'}, RES_All.Name);
clear drugs i

% ---- Population conditions & drug mapping ------------------------------
puc_base = @(puc) true(height(puc),1);

PopulationConditions = cell2table({
    "base"  "RF"   "isGoodRF"       puc_base "all" [0 365]
    "F5080" "RF"   "isGoodRF_F5080" puc_base "all" [0 365]
    "base"  "PHEN" "isPHEN"         puc_base "all" [0 365]
    }, 'VariableNames', {'tag','type','pop_mask','puc_mask','puc_time_method','timeframe'});
nCond = height(PopulationConditions);

DrugMapping = cell2table({ ...
    "RF"   "\beta-lactam"  "cov_Drug"       "beta_lactam"     "beta_lactam"        "cov"
    "RF"   "FOS"           "cov_Drug"       "fosfomycin"      "fosfomycin"         "cov"
    "RF"   "TMP+SMX"       "cov_TMP_SMX"    ""                "Trimethoprim_Sulfa" "cov"
    "RF"   "CPR (gene)"    "cov_Drug"       "quinolone"       "quinolone"          "cov"
    "RF"   "CPR (mut)"     "mut_Drug"       ""                "quinolone"          "mut"
    "PHEN" "CEF"           "clean_counts_n" "CEF_A"           "CEF"                "phen"
    "PHEN" "NIT"           "clean_counts_n" "NIT-8_B"         "Nitrofurantoin"     "phen"
    "PHEN" "FOS"           "clean_counts_n" "FOS-16_B"        "fosfomycin"         "phen"
    "PHEN" "TMP"           "clean_counts_n" "TMP_A"           "Trimethoprim_Sulfa" "phen"
    "PHEN" "CPR"           "clean_counts_n" "CPR_A"           "quinolone"          "phen"
    }, 'VariableNames', {'type','pretty','entity','field','ab_group','res_type'});

DayEdges = [-inf, -1.5:365.5, inf];

% ---- Build patient-paired PUC x FOBT data ------------------------------
CondResults = struct;
for iCond = 1:nCond
    cond_row = PopulationConditions(iCond,:);
    cr = build_matching_fobt_puc_data(DrugMapping, cond_row, ...
        POP, Cohorts, PUC, MHS_ResGroups, I_as_RS, DayEdges, RF, RF_Entities, PHEN, PHEN_Entities);
    CondResults.(cr.label) = cr;
end

res_th_type2params = cell2table({...
    'cov'  Res_th.cov  [1e-4 1] 1
    'mut'  Res_th.mut  [0 4]    0
    'phen' Res_th.phen [1e-2 1] 1
    },'VariableNames',{'res_type','res_th','x_lim','islog'}, 'RowNames',{'cov','mut','phen'});

% ---- CDF of faecal resistance stratified by urine-culture resistance ---
cond_labels = string(fields(CondResults));
for iCond = 1:numel(cond_labels)
    cond_label = cond_labels(iCond);
    cr = CondResults.(cond_label);
    ConfMatDrugs = nan(2,2,cr.nDrugs);
    fprintf('\n%25s %s\n', '"'+cond_label+'"', repmat('-',1,50));
    fobt_ths = nan(1, cr.nDrugs);
    for ijD = 1:cr.nDrugs
        jDrug = cr.jDrugs(ijD);
        dmap = table2struct(DrugMapping(jDrug,:));
        param = res_th_type2params(dmap.res_type,:);
        if dmap.pretty == "\beta-lactam", nres = 3; else, nres = 1; end
        fobt_th = param.res_th;
        fobt_ths(ijD) = fobt_th;

        [~, num_measured, num_res] = count_UTI_resistances(...
            cr.v123count(:,ijD,:,:), DayEdges, cr.timeframe, cr.puc_time_method);
        fobt_val = cr.fobt_val(:,ijD);

        [conf_mats, chosen, Fisher_p, KS_p] = ...
            get_fobt_puc_confusion_mat_and_association_statistics(...
            fobt_val, fobt_th, num_measured, num_res, nres);
        fprintf('%25s  Fisher=%s\n', dmap.pretty, sprintf('%6.4f ', Fisher_p));
        ConfMatDrugs(:,:,ijD) = conf_mats(:,:,1);

        if any(sum(chosen,1) < 4)
            continue
        end
        if any(cond_label == ["RF_base" "PHEN_base"])
            fig_num = 13000 + iCond*100 + ijD;
            [cf, h] = fig_cdf_res_sen(fig_num, cr.type, fobt_val, chosen, dmap.pretty, ...
                fobt_th, param.islog, param.x_lim, CLRS.sen, CLRS.res);
            if ~isempty(h)
                print_figure(cf, sprintf('fig_cdf_res_sen_%s_%s', cond_label, strrep(dmap.pretty,'\','')));
            end
        end
    end
    CondResults.(cond_label).fobt_ths = fobt_ths;
    CondResults.(cond_label).ConfMatDrugs = ConfMatDrugs;
end

% ---- Risk bar charts ----------------------------------------------------
for cond_label = ["RF_base" "PHEN_base"]
    cr = CondResults.(cond_label);
    rf_or_phen = cr.type;
    drugs = cellstr(DrugMapping.pretty(DrugMapping.type == rf_or_phen));
    conf_mat = cr.ConfMatDrugs;                 % [2,2,nDrugs]
    conf_mat = permute(conf_mat, [1 2 4 3]);    % [2,2,1,nDrugs]
    [cf, h] = fig_bar_res_sen(10200, conf_mat, {{}, {}, drugs}, rf_or_phen, CLRS.sen, CLRS.res, []);
    if ~isempty(h)
        print_figure(cf, sprintf('fig_bar_res_sen_%s', cond_label));
    end
end

% ---- Odds ratio vs time (wrapped: heaviest / most data-hungry step) ----
try
    rng(54)
    for cond_label = "RF_base"
        cr = CondResults.(cond_label);
        jd = find(DrugMapping.pretty(cr.jDrugs) ~= "CPR (gene)");
        F3 = convert_fobt_val_to_F3(cr.fobt_val, cr.fobt_ths);
        F3 = F3(:,jd);
        v123 = cr.v123count(:,jd,:,:);
        assert(DayEdges(1)==-inf & DayEdges(end)==inf)
        day_edges = DayEdges(2:end-1);
        v123 = v123(:,:,2:end-1,:);
        days = (day_edges(1:end-1) + day_edges(2:end)) / 2;
        V4 = convert_v123count_to_v4(v123, [], cr.puc_time_method, 3);
        jT = find(days==cr.timeframe(1)):56:365;
        days = days(jT);
        V4 = V4(:,:,jT);
        [~, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile] = ...
            calc_F3_V4_association_across_time_and_drugs(F3, V4, days, [16 50 84], 'complete', 2000, 2000);
        cf = fig_oddsratio_vs_time(30, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile, days);
        print_figure(cf, sprintf('fig_oddsratio_vs_time_%s', cond_label));
    end
catch ME
    fprintf('WARNING: odds-ratio figure skipped (%s)\n', ME.message);
end

% ---- Resistance-gene coverage heatmap ----------------------------------
ok = RF.isGood & RF.FOBTisFirst & ~RF.isrep;
[cf, gene_cov_map.num_samples_not_plotted] = ...
    fig_res_gene_coverage_matrix(101, RF(ok,:), RF_Entities, 'Gene', 'Drug', Res_th.cov);
record_value('Gene coverage heatmap', [], gene_cov_map)
print_figure(cf, 'fig_res_gene_coverage_matrix', [], [], 'Arial')

% ---- Population statistics ---------------------------------------------
ValidSMP = SMP(inrange(SMP.DateDiff, Timespan.BAC_FOBT_inclusion) & ...
    ismember(SMP.FOBT_ID, POP.FOBT_ID), :);
cf = fig_pop_stat(102, POP, ValidSMP, Cohorts.Coinciding, Cohorts.NonCoinciding);
print_figure(cf, 'fig_pop_stat');

% ---- Write recorded values ---------------------------------------------
values = record_value();
writetable(values, 'values.xls');

fprintf('\nDemo complete. Outputs in: %s\n', out_dir);
fprintf('  figures/  (SVG files)\n  values.xls\n');

end
