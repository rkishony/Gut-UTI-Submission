function varargout = singleton_expand(varargin)
% Expand inputs along singleton dims to a common size

% Get common tgt size
maxN  = max(cellfun(@ndims, varargin));
sizes = cellfun(@(x) size(x,1:maxN), varargin, 'uni', false);

tgt   = max(cat(1, sizes{:}), [], 1);

% Verify expandability:
for k = 1:numel(sizes)
    s = sizes{k};
    if any((s ~= 1) & (s ~= tgt))
        error('Inputs can only be expanded along singleton dimensions.');
    end
end

% Singleton expansion:
varargout = cellfun(@(x,s) repmat(x, tgt./s), varargin, sizes, 'uni', false);

end
