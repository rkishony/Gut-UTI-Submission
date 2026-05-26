function varargout = load_or_create(scriptname_or_func, filename, force_create, varargin)
% Loads MAT file. If missing, runs script to create it.
% If no output is requested, assigns vars to caller like `load`.
if nargin<3
    force_create = false;
end

[~,dispfile, ext] = fileparts(filename);
dispfile = [dispfile ext];

f = which(filename);
[date, date_str] = get_date(f);
if isempty(f)
    op = 'CREATE';
elseif force_create
    op = 'RE-CREATE';
else
    op = 'LOAD';
end

is_func = isa(scriptname_or_func, 'function_handle');
if is_func
    func_name = func2str(scriptname_or_func);
else
    func_name = scriptname_or_func;
end

if strcmp(op, 'LOAD')
    fprintf('File %26s exists (%s), Loading ...\n', add_quotes(dispfile), date_str)
else
    fprintf('======   running %26s to %10s %26s ======\n', add_quotes(func_name), op, add_quotes(dispfile))
    tic
    if is_func
        [varargout{1:nargout}] = scriptname_or_func(varargin{:});
        cachedData = varargout;
        save(filename,'cachedData');
    else
        out_filepath = run_script(scriptname_or_func);

        date2 = get_date(f);

        if isnan(date2) || date==date2
            error('File was not re-created')
        end
    end
    fprintf('====== completed %26s in %6.1f sec %26s ======\n\n', add_quotes(func_name), toc, '')
end

L = load(filename);

if is_func
    varargout = L.cachedData;
else
    if nargout == 0
        vars = fieldnames(L);
        for i = 1:numel(vars)
            assignin('caller', vars{i}, L.(vars{i}));
        end
        varargout = {};
    else
        varargout = {L};
    end
end
end

function out_filepath = run_script(scriptname)
% we assume that out_filepath is created by the script.
run(scriptname);
end
      
function [date, date_str] = get_date(filepath)
if exist(filepath, 'file')
    d = dir(filepath);
    date = d.datenum;
    date_str = datestr(date);
else
    date = nan;
    date_str = 'MISSING';
end
end

function name = add_quotes(name)
name = ['`' name '`'];
end