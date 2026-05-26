function folder = get_folders(name)

base = fileparts(fileparts(fileparts(mfilename('fullpath'))));

folders.analyses_server = fullfile(base, 'output');
folders.figures = fullfile(base, 'figures');


if nargin == 0
    folder = folders;
else
    folder = folders.(name);
end

end
