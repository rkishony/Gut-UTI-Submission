function add_required_paths()

here = fileparts(mfilename('fullpath'));

cwd = pwd;
cd(fullfile(here, 'analyses'))
folders = get_folders();
cd(cwd)

addpath(genpath(folders.loading_scripts));
addpath(genpath(folders.utils));
addpath(genpath(folders.analyses));

end