function [PHEN, PHEN_Abs] = load_CommunityResPhen(folder_name)


note_mapping = {
    '',          0
    '+haze',     inf
    '+low haze', 0
    'haze',      inf
    'lawn',      inf
    'lawn?',     inf
    'low haze',  0
    'merged',    nan  
    'mid haze',  0
    };
  

%% Set out_filepath:
folders = get_folders();

% Read in table as strings
filename_comm_exp = fullfile(folders.from_Dropbox, folder_name, 'results.csv');
% opts = detectImportOptions(filename_comm_exp);
% opts = setvartype(opts,opts.VariableNames,"string");
PHEN = readtable(filename_comm_exp,'VariableNamingRule','preserve');
PHEN(:,{'well', 'plate_num', 'isolate_id'}) = [];

vars = string(PHEN.Properties.VariableNames);
PHEN.Properties.VariableNames{vars == "name"} = 'Sample';
PHEN.Properties.VariableNames{vars == "name"} = 'Sample';

count_prefix = 'AllAbs_NumColonies_';
is_count = startsWith(vars, count_prefix);
PHEN_Abs = cellfun(@(x) x(numel(count_prefix)+1:end),vars(is_count),'UniformOutput',false);
counts = PHEN{:, is_count};

note_prefix  = 'AllAbs_Note_';
is_note = startsWith(vars, note_prefix);

notes = PHEN{:, is_note};
notes = lower(notes);
unique_notes = unique(notes);

assert(isempty(setdiff(unique_notes, note_mapping(:,1))))

%%
% Process the measurement string to define 'counts', 'ismerged', 'islawn'
for j = 1:size(note_mapping,1)
    [key, val] = note_mapping{j, :};
    match = strcmp(notes, key);
    counts(match) = counts(match) + val;
end

PHEN.counts = counts;
PHEN.ismerged = strcmp(notes, 'merged');

% Process the Sample name to define 'isU', 'isrep', 'isempty'
PHEN.isUTI = PHEN.UTI_ID > 0;
PHEN.isempty = strcmp(PHEN.Sample,'Neg. Ctrl');
PHEN.isFOBT = ~(PHEN.isempty | PHEN.isUTI);
PHEN = PHEN(PHEN.isFOBT,:);

% Get sample FOBT_ID and UTI_ID
PHEN = sortrows(PHEN, {'isFOBT', 'FOBT_ID', 'rep'}, {'descend', 'ascend', 'ascend'});

% Assert only one FOBT_ID per ~isrep
is_first = [true; PHEN.FOBT_ID(2:end)~=PHEN.FOBT_ID(1:end-1)];
wrong_rep = is_first & PHEN.rep & ~isnan(PHEN.FOBT_ID);
if any(wrong_rep)
    fprintf('WARNING: %d PHEN are mistakenly labeled as isrep. Correcting\n', sum(wrong_rep))
    PHEN.isrep(wrong_rep) = false;
end
wrong_nonrep = sum(~is_first & ~PHEN.rep & ~isnan(PHEN.FOBT_ID));
if any(wrong_nonrep)
    fprintf('WARNING: %d PHEN are mistakenly labeled as ~isrep. Correcting\n', sum(wrong_nonrep))
    PHEN.isrep(wrong_nonrep) = true;
end

end
