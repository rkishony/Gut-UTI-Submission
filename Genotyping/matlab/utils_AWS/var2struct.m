function s = var2struct(varargin)
for k = 1:nargin
    nm = inputname(k);
    if isempty(nm)
        error('Argument %d must be a named variable.', k);
    end
    s.(nm) = varargin{k};
end
end
