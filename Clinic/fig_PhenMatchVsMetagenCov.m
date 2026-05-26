
% THIS SCRIPT takes as input
% - isolate phenotyping data
% - coverage breadth data of UTI isolates by pathobiome metagenomic reads
% - Maccabi urine cultures data

% A datapoint is a UTI isolate with a coverage breadth associated with it,
% we color the datapoint by whether or not the UTI matches the FOBT
% phenotypically, where a match is calculated according to the minimal
% distance between the UTI (MHS-measured profile) and FOBT isolates from the 
% same patients. 
% The FOBT data is binary and pooled from three experiments.


add_required_paths
folders = get_folders;
addpath(genpath(folders.IsolatePhenotyping))
folders = get_folders;
if ~exist('CONGIF','var')
    global CONFIG
    CONFIG.print_figures = true;
end
record_value([]); % clear result recording table

calculate_Pvalue = true;

% load phenotypic data (From dropbox):
all_res = readtable(fullfile(folders.from_Dropbox, 'IsolateResData.csv'));

% load genotypic data
alignment_table = load_alignment;

% etha0 is the fraction uncovered -- uncov_len./genome_len.

%%

% choose a range of number of isolates per patient. [0 100] includes any number of isolates.
r = [0 100]; 


%% process genotypic data

etha0_table = alignment_table(:,{'uti_sample_name','UTI_ID','etha0'}); 

% take only isolate A, parse out patient ID
for sample = 1:height(etha0_table)
    sample_name = etha0_table.uti_sample_name{sample};
    parts = strsplit(sample_name,'.');
    etha0_table.patient_ID(sample) = str2double(parts{1}(2:end));
    etha0_table.isolate{sample} = parts{3};
end

unique_etha0_table = table;
count = 0
considered_isolate = false([1, height(etha0_table)]);
for u = unique(etha0_table.UTI_ID)'
    count = count+1;
    ind = find(strcmp(etha0_table.isolate,'A') & etha0_table.UTI_ID==u);
    if isempty(ind)
        ind = find(strcmp(etha0_table.isolate,'B') & etha0_table.UTI_ID==u);
    end

    considered_isolate(ind) = true;
end

etha0_table = etha0_table(considered_isolate,:);

% assert that all UTIs are present and only once
assert(height(etha0_table)==length(unique(alignment_table.UTI_ID)))

% table --> vectors
etha0 = etha0_table.etha0;
uti_gid = etha0_table.UTI_ID;
fobt_gid = etha0_table.patient_ID;

record_value('Phenotyping analysis cohorts', 'ana_cohorts.MatchPhenGen.initial_uti_cultures', length(unique(uti_gid)))
record_value('Phenotyping analysis cohorts', 'ana_cohorts.MatchPhenGen.initial_uti_patients', length(unique(fobt_gid)))

%% load/process Maccabi data
redo = true;
if ~exist(fullfile(folders.mat_outputs, 'BACres_complexUTIssplit.mat'),'file') | redo
    
    load(fullfile(folders.mat_outputs, 'allMaccabiData.mat')); 

    MHSall.BAC = BAC;
    MHSall.RES = RES_All;
    MHSall.UTIpop = UTIpop;

    MHSdata = processMHSdata(MHSall, unique(uti_gid));
else
    load(fullfile(folders.mat_outputs, 'BACres_complexUTIssplit.mat'));
end

%% process phenotyping data:

ABs = all_res.Properties.VariableNames(end-4:end); % 5 antibiotics, the last 5 columns in the table

resistancemaps = table2array(all_res(:,ABs));

%% Which patients have both a UTI isolate and an FIT isolate in the data?

valid = all_res.rep==0 & ~all_res.from_lawn & ~all_res.is_neg_ctrl & all_res.FOBT_ID>0;
patients_w_fobt = unique(all_res.FOBT_ID(all_res.FOBT_ID>0 & all_res.UTI_ID==0 & valid));


% CHOOSE:
patients = patients_w_fobt; 

npatients = length(patients); %length(pairedpatients);

%% Create the patients' personal matrices of resistance:
personalFOBT_inds = cell(npatients, 1);

personalmats = cell(npatients, 1);

for p = 1:npatients
    personalFOBT_inds{p} = find(patients(p)==all_res.FOBT_ID & all_res.UTI_ID==0 & valid);

    personalmats{p} = resistancemaps(personalFOBT_inds{p},:);

end

%% analyse
rel_patients = intersect(patients, fobt_gid);
rel_personalmats = personalmats(ismember(patients,rel_patients));

etha0s_and_distances_data = analyze_relation_between_UTIcovbre_and_FITprofiles([etha0, uti_gid, fobt_gid], MHSdata, rel_personalmats, rel_patients);
[hbin, binedges] = breadths_and_distances_into_colorcoded_bars_data(etha0s_and_distances_data, r);

record_value('Phenotyping analysis cohorts', 'ana_cohorts.MatchPhenGen.valid_uti_patients', length(unique(etha0s_and_distances_data.patientID)))
record_value('Phenotyping analysis cohorts', 'ana_cohorts.MatchPhenGen.valid_uti_cultures', sum(hbin,'all'))

identical_rate_at_genotypic_match = hbin(1,1)/sum(hbin(1,:));
%% plot

plot_IsolatePhenMatch_vs_CovBreadth(hbin, binedges, 'PhenMatchV_vs_MetagenCov')

print_figure(gcf, 'PhenMatchV_vs_MetagenCov')

%% pval through permutation test
if ~calculate_Pvalue
    return
end

% we permute the breadths and seehow often we get a proportion of
% identicals equal or higher in the fully-covered column
nperm = 10000;
for n = 1:nperm
    permuted_breadths_and_distances_data = etha0s_and_distances_data;
    perm = randperm(height(etha0s_and_distances_data));
    permuted_breadths_and_distances_data.etha0 = etha0s_and_distances_data.etha0(perm);
    permuted_breadths_and_distances_data.capped_etha0 = etha0s_and_distances_data.capped_etha0(perm);

    [permuted_hbin(:,:,n), permuted_binedges(:,:,n)] = breadths_and_distances_into_colorcoded_bars_data(permuted_breadths_and_distances_data, r);
end

pval = mean(permuted_hbin(1,1,:)>=hbin(1,1));
record_value('pval breadth is associated with phen profile match', [], pval)


values = record_value();
writetable(values, 'phen_values.xls', 'WriteMode','replacefile')
