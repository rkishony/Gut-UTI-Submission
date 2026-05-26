function T = loadisolateMLST()

folders = get_folders();
summary_file = fullfile(folders.analyses_server, 'summary', 'mlst_isolates.csv');

T = readtable(summary_file, 'TextType', 'string');

% logicals
T.is_faecal = strcmpi(string(T.is_faecal), "True");
T.is_isolate = strcmpi(string(T.is_isolate), "True");

% normalize species names
num_species = numel(unique(T.species));
T.species = normalize_species(T.species);
assert(num_species == numel(unique(T.species)));

% assert all are isolates
assert(all(T.is_isolate));
T.is_isolate = [];

end


function name = normalize_species(name)
% remove everything after the first '_' for each entry
name = string(name);
has_sep = contains(name, "_");
if any(has_sep)
    parts = split(name(has_sep), "_");
    name(has_sep) = parts(:, 1);
end
end
