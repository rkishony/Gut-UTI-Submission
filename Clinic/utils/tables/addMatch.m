function T = addMatch(T, S, key, fields, new_names)
if nargin < 5, new_names = []; end
if ischar(fields), fields = cellstr(fields); end

[~, idx] = unique(S(:, key), 'stable');

sub = S(idx, [key, fields]);

% rename new_names is provided
if ~isempty(new_names)
    if ischar(new_names), new_names = cellstr(new_names); end
    sub.Properties.VariableNames(2:end) = new_names;
end
T.indx = (1:height(T))';
T = outerjoin(T, sub, 'Keys', key, 'Type', 'left', 'MergeKeys', true);
T = sortrows(T, 'indx');
T.indx = [];
end
