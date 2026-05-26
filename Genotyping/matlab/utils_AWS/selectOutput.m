function out = selectOutput(func, idx, varargin)
outs = cell(1, idx);
[outs{:}] = func(varargin{:});
out = outs{end};
end