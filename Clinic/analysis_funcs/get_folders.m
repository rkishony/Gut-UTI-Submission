function folders = get_folders()

[current, ~, ~] = fileparts(mfilename('fullpath'));
root = fileparts(current);

folders.project_root = root;
folders.loading_scripts = fullfile(root, 'loading_scripts');
folders.fig_scripts = fullfile(root, 'fig_scripts');
folders.translation_tables= fullfile(root, 'source_data', 'translation_tables');
folders.utils = fullfile(root, 'utils');
folders.analysis_funcs = fullfile(root, 'analysis_funcs');

folders.Maccabi_data = fullfile(root, 'source_data', 'Maccabi_data');
folders.collection = fullfile(root, 'source_data', 'uti_fobt_collection');
folders.mat_outputs = fullfile(root, 'mat_outputs');
folders.meta_data= fullfile(root, 'source_data', 'meta_data');
folders.from_Dropbox = fullfile(root, 'source_data', 'from_Dropbox');
folders.from_labserver = fullfile(root, 'source_data', 'from_labserver');
folders.IsolatePhenotyping = fullfile(root, 'IsolatePhenotyping');

folders = structfun(@resolvePath, folders, 'UniformOutput',false);
end


function absPath = resolvePath(relPath)
    absPath = char(java.io.File(relPath).getCanonicalPath());
end