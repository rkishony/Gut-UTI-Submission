function f = get_folders(folder_name)

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
assert(endsWith(root, 'Manuscript-phen-only'))

f = struct;
f.root = root;
f.utils = fullfile(root, 'utils');
f.loading_scripts = fullfile(root, 'loading_scripts');
f.analyses = fullfile(root, 'analyses');
f.analyses_server = fullfile(root, 'analyses_server');
% f.mat_outputs = fullfile(root, 'mat_outputs');
f.metadata = fullfile(root, 'metadata');
f.source_data = fullfile(root, 'source_data');
f.figures = fullfile(root, 'figures');
f.image_analysis = fullfile(root, 'PHEN_raw_and_preprocessing', 'image_analysis');
f.raw_images = fullfile(root, 'PHEN_raw_and_preprocessing', 'raw_data');
f.from_MHS = fullfile(root, 'from_MHS');
f = structfun(@resolvePath, f, 'UniformOutput',false);

mkdir_if_needed(f.figures)
% mkdir_if_needed(f.mat_outputs)
mkdir_if_needed(f.source_data)

if nargin>0
    f = f.(folder_name);
end



end

function mkdir_if_needed(pth)
if ~exist(pth, 'dir')
    mkdir(pth)
end
end
function absPath = resolvePath(relPath)
    absPath = char(java.io.File(relPath).getCanonicalPath());
end
