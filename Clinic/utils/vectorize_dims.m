function Y = vectorize_dims(fun, X, dims)
%VECTORIZE_DIMS   Apply FUN to sub-arrays of X spanning dimensions DIMS.
otherDims = setdiff(1:ndims(X), dims);
C         = num2cell(X, otherDims);            % cell array of size X(dims)
cellSz    = size(C);
Yc        = cellfun(fun, C, 'UniformOutput', false);
outSz     = size(Yc{1});
Yall      = cat(numel(outSz)+1, Yc{:});        % stack along new dim
Y         = reshape(Yall, [outSz, cellSz]);
end
