function [M, xs, ys] = table_pivot(func, default, tbl, zcol, xs, ys)
% PIVOT  Create M(i,j)=func(z) for each unique x,y.

xcol = xs.Properties.VariableNames;
ycol = ys.Properties.VariableNames;

xv = tbl(:, xcol);
yv = tbl(:, ycol);
zv = tbl.(zcol);

[~, ix] = ismember(xv, xs, 'rows');
[~, iy] = ismember(yv, ys, 'rows');

M = accumarray([ix, iy], zv, [size(xs,1), size(ys,1)], func, default);
end
