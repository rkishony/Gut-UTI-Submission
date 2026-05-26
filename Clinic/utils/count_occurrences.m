function counts = count_occurrences(x, varargin)
[~, ~, ic] = unique(x, varargin{:});
freq = accumarray(ic, 1);
counts = freq(ic);
end
