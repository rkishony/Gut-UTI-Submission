function T = combtable(S)
%COMBTABLE All combinations of column values in struct S -> table.
%   Fields of S are variable names; each field is a vector/cell/strings.

names = fieldnames(S);
vals  = cellfun(@(f) S.(f)(:), names, 'uni', false);
n     = cellfun(@numel, vals);
idx   = arrayfun(@(m) 1:m, n, 'uni', false);
I     = cell(1,numel(idx));
[I{:}] = ndgrid(idx{:});
cols  = cellfun(@(v,Ik) v(Ik(:)), vals', I, 'uni', false);
T     = table(cols{:}, 'VariableNames', names);
end
