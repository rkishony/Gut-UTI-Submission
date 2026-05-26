cd analysis_funcs/
folders = get_folders();
cd ..

addpath(genpath(folders.utils))
addpath(folders.analysis_funcs)
addpath(folders.loading_scripts)
addpath(folders.fig_scripts)

addpath(folders.mat_outputs);
addpath(genpath(folders.IsolatePhenotyping))
