function add_required_paths()

base = fileparts(mfilename('fullpath'));

addpath(fullfile(base, 'loading_scripts'));
addpath(genpath(fullfile(base, 'utils_AWS')));
addpath(genpath(fullfile(base, 'utils')));
addpath(genpath(fullfile(base, 'figures')));

end
