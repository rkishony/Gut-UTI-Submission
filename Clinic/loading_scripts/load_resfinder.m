% auto-script. Runs by `load_or_create`. Can also run manually from the command line.

%% Set out_filepath:
clear
folders = get_folders();
out_filepath = fullfile(folders.mat_outputs, 'resfinder.mat');


%% folders and sample names
server_folder = folders.from_labserver;
resfinder_dir = fullfile(server_folder, 'resfinder');
TRIM = load_trimstats();

RESFDR = table;
RESFDR.Sample = TRIM.run_name;
RESFDR.isFOBT = TRIM.is_faecal;
RESFDR.FOBT_ID = TRIM.fid;
RESFDR.isrep = TRIM.is_rep;
is_downsampled = true;  % For ResFinder we used the full data, not the downsampled.
if is_downsampled
    RESFDR.ReadsTotal = TRIM.effective_num_reads;  
    RESFDR.bpTotal = TRIM.effective_total_bases;
else
    RESFDR.ReadsTotal = TRIM.num_reads;  
    RESFDR.bpTotal = TRIM.total_bases;
end
RESFDR = RESFDR(TRIM.is_faecal & ~TRIM.is_isolate,:);

%% 

RESFDR = sortrows(RESFDR, {'isFOBT', 'FOBT_ID', 'isrep'}, {'descend', 'ascend', 'ascend'});
RESFDR.jRF = (1:height(RESFDR))';

%% get gene coverage

% define drugs
drugfiles = dir([resfinder_dir filesep '*' filesep 'acq' filesep 'resfinder_kma' filesep 'kma_*.res']);
udrugfiles = unique({drugfiles.name});
drugs = cellfun(@(x) x(strfind(x,'_')+1:strfind(x,'.res')-1),udrugfiles,'UniformOutput',false);

% Standardize drug names to allow using as matlab stuct fields:
new_drug_names = cellfun(@makeValidName, drugs, 'uniformoutput', false);  % 'beta-lactam' -> 'beta_lactam'

% read .res files
resfdr_cov_names = {'Sample','Drug','Template','Score','Expected','Template_length','Template_Identity',...
    'Template_Coverage', 'Query_Identity', 'Query_Coverage', 'Depth', 'q_value', 'p_value'};
RESFDR_COV = table;

fprintf('%45s', 'Getting resistance-gene coverage per sample ')
n_samples = height(RESFDR);
for isample = 1:n_samples
    if ~mod(isample,10), fprintf('.'); end
    for idrug = 1:length(drugs)
        filename = fullfile(resfinder_dir, RESFDR.Sample{isample}, 'acq', 'resfinder_kma', ['kma_' drugs{idrug} '.res']);
        if exist(filename,'file')
            rf_cov = readtable(filename,'FileType','text',...
                'ReadVariableNames',false,'HeaderLines',1,'Delimiter','tab');
            h = height(rf_cov); 
            if h > 0
                rf_cov.Properties.VariableNames = resfdr_cov_names(3:end);
                rf_cov{:,'Sample'} = RESFDR.Sample(isample);
                rf_cov{:,'Drug'} = new_drug_names(idrug);
                rf_cov = rf_cov(:,resfdr_cov_names);
                rf_cov.Singleton = repmat(h==1, h, 1);
                RESFDR_COV = [RESFDR_COV; rf_cov];
            end
        end
    end
end
fprintf('\n')

%% read all PointFinder_results.txt
pointfvarnames = {'Sample','Gene','AAlocus','Resistance','PMID'};
RESFDR_MUT = table;
fprintf('%45s', 'Getting resistance-mutations per sample ')
for isample = 1:n_samples
    if ~mod(isample,10), fprintf('.'); end
    filename = fullfile(resfinder_dir, RESFDR.Sample{isample}, 'pointEcoli', 'PointFinder_results.txt');
    if exist(filename,'file')
        warning off
        rf_mut = readtable(filename,'FileType','text',...
            'ReadVariableNames',true,'HeaderLines',0,'Delimiter','tab');
        warning on
        h = height(rf_mut);
        if h > 0
            rf_mut{:,'Sample'} = RESFDR.Sample(isample);
            rf_mut.Singleton = repmat(h==1, h, 1);
            RESFDR_MUT = [RESFDR_MUT; rf_mut];
        end
    end
end
fprintf('\n')
RESFDR_MUT(:,'Gene') = cellfun(@(x) x(1:strfind(x,'.')-3),RESFDR_MUT.Mutation,'UniformOutput',false);
RESFDR_MUT(:,'AAlocus') = cellfun(@(x) x(strfind(x,'.')-1:end),RESFDR_MUT.Mutation,'UniformOutput',false);
RESFDR_MUT = RESFDR_MUT(:,[pointfvarnames, 'Singleton']);

% Fix resistance names:
RESFDR_MUT.Resistance(strcmp(RESFDR_MUT.Resistance, 'Nalidixic acid,Nalidixic acid,Ciprofloxacin')) = {'Nalidixic acid,Ciprofloxacin'};

%% Standardize field naming of RF_MUT and RF_COV:
RESFDR_COV = renamevars(RESFDR_COV, 'Template', 'Allele'); 
RESFDR_COV.Gene = cellfun(@(x) regexp(x, '^[a-zA-Z]+','match','once'), RESFDR_COV.Allele, 'UniformOutput',false);

RESFDR_MUT = renamevars(RESFDR_MUT, {'Resistance', 'AAlocus'}, {'Drug', 'Allele'});
RESFDR_MUT.Depth = ones(height(RESFDR_MUT),1);  % to be able to count similarly to RF_COV


%% save resTable
save(out_filepath, 'RESFDR', 'RESFDR_COV', 'RESFDR_MUT')
