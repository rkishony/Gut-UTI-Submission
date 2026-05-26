function funcMem = file_memorize(func, force_create)
if nargin<2
    force_create = false;
end
info = functions(func);
filename = info.function;

folders = get_folders();
filepath = fullfile(folders.mat_outputs, [filename '.mat']);

funcMem = @wrapper;

function varargout = wrapper(varargin)
[varargout{1:nargout}] = load_or_create(func, filepath, force_create, varargin{:});
end
end