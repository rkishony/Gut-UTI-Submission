function tbl = extractCols(tbl, cols, newNames)
if isempty(cols)
    cols = 1:size(tbl, 2);
end
tbl = tbl(:, cols);
if nargin>=3
    tbl.Properties.VariableNames = newNames;
end
end