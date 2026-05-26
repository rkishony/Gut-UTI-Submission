%% Setup
close all
clear all
clc
add_required_paths()
record_value([]);  % clear result recording table

RERUN_ALL_IMAGE_ANALYSIS = false;

folders = get_folders();


%% Analyse all phenotype images and create PDF and XLS outputs
tic
exp_folders = {'CommunityPhenotyping', 'CommunityPhenotypingB','IsolateResA', 'IsolateResB', 'IsolateResC'};
base_folder = get_folders('image_analysis');

destination_folder = fullfile(get_folders('source_data'), 'Phenotyping');
dlt_folder_contents(destination_folder)
for iExp = 1:numel(exp_folders)
    exp_folder = exp_folders{iExp};
    analysis_path = fullfile(base_folder, exp_folder);
    destination_folder_specific = fullfile(destination_folder, exp_folder);
    mkdir(destination_folder_specific)

    if RERUN_ALL_IMAGE_ANALYSIS
        filepath = fullfile(analysis_path, 'analyzeRES_here');
        runScriptAsFunc(filepath)
    end

    copyfile(fullfile(analysis_path, 'PlateLayout.xlsx'), destination_folder_specific)
    if exist(fullfile(analysis_path, 'input', 'readme.xlsx'),'file')
        copyfile(fullfile(analysis_path, 'input', 'readme.xlsx'), destination_folder_specific)
    end
    copyfile(fullfile(analysis_path, 'input', 'strips.xlsx'), destination_folder_specific)

    copyfile(fullfile(analysis_path, 'output', '*.*'), destination_folder_specific)
end
toc


%% Combine XLS outout of phen resistance from Exp A, B, C.

all_res = combineExperiments();


%% Calculate phenotypic distance statisitcs

distinguishing_conc = readtable(fullfile(get_folders('metadata'), 'distinguishing_conc.xlsx'));
Abs = distinguishing_conc.Properties.VariableNames;

[b_same, b_bs, nprofs, nperm, by_chance_in_nperm] = phenotypic_distance_analysis(all_res, Abs);
assert(all(nprofs>0))

pval_permutations_with_match_ratio_same_or_higher = by_chance_in_nperm/nperm;
if pval_permutations_with_match_ratio_same_or_higher==0
    pval_permutations_thresh = (by_chance_in_nperm+1)/nperm;
    record_value('Phenotyping analysis', 'pval_permutations_thresh', pval_permutations_thresh)
end

record_value('Phenotyping analysis', 'pval_permutations_with_match_ratio_same_or_higher', pval_permutations_with_match_ratio_same_or_higher)
record_value('Phenotyping analysis', 'total_permutations_with_match_ratio_same_or_higher', nperm)

%% Figure of phenotypic profile match self versus other

fig_hist_phen_profiles_match_same_vs_diff_pat(1011, b_same, b_bs, 32);
print_figure(gcf, fullfile(get_folders('figures'), ...
    'fig_hist_phen_profiles_match_same_vs_diff_pat'))


%% Figure of diversity of phenotypic profiles per patient

h = fig_hist_num_phen_res_fobt_iso_profiles_per_pat(1012, nprofs);

print_figure(gcf, fullfile(get_folders('figures'), ...
    'fig_hist_num_phen_res_fobt_iso_profiles_per_pat'))

MultiProfPatientPerc = 100*mean(nprofs>1);
record_value('Phenotyping analysis', 'Perc of patients with >1 resistance profile', MultiProfPatientPerc)
MeanNumResProfs = mean(nprofs);
record_value('Phenotyping analysis', 'Mean number of resistance profiles per patient', MeanNumResProfs)

%% Calculate and record patient cohort stat

filepath = fullfile(folders.analyses, 'IsolatePhenotyping', 'output_patient_and_isolate_sets');
[xp_cohorts, analysis_cohorts] = runScriptAsFunc(filepath,{'S','A'});

exp_cohorts = struct;
flds = fields(xp_cohorts);
for i = 1:numel(flds)
    fld = flds{i};
    if numel(xp_cohorts.(fld))==1 && isa(xp_cohorts.(fld),'double')
        exp_cohorts.(fld) = xp_cohorts.(fld);
    end
end
record_value('Phenotyping experiment cohorts', [], exp_cohorts)

ana_cohorts = struct;
flds = fields(analysis_cohorts);
for i = 1:numel(flds)
    fld = flds{i};
    ana_cohorts.(fld) = analysis_cohorts.(fld).total;
end
record_value('Phenotyping analysis cohorts', [], ana_cohorts)

% record: number of isolates excluded for no-growth in mc:
isMeasured = is_measured(all_res);
valid = is_valid(all_res, isMeasured);
n_invalid_utis = height(all_res(isMeasured & ~valid & all_res.UTI_ID>0,:));
n_invalid_fits = height(all_res(isMeasured & ~valid & all_res.UTI_ID==0,:));

record_value('Phenotyping analysis', 'Number of UTI isolates with no-growth in MC', n_invalid_utis) % note this is isolates like in the paper, not cultures like the ones recorded in summarize_isolates_in_experiments.m
record_value('Phenotyping analysis', 'Number of FIT isolates with no-growth in MC', n_invalid_fits)

%% Save all recorded values

values = record_value();
writetable(values, 'phen_values.csv')
