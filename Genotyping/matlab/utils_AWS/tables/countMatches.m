function [countsA, countsB] = countMatches(A, B, keyVars, condFcn, bool)
% Counts matching rows between two tables
if nargin<4
    condFcn = [];
end

if nargin<5
    bool = false;
end

nA = height(A);
if bool
    countsA = zeros(nA,1,'logical');
else
    countsA = zeros(nA,1);
end

for i = 1:nA
    maskKey = ismember( B(:, keyVars), A(i, keyVars), 'rows' );
    if ~isempty(condFcn)
        Bsub    = B(maskKey, :);
        maskKey = condFcn( A(i,:), Bsub );
    end
    countsA(i) = sum(maskKey);
end

if nargout<2
    return
end

nB = height(B);
if bool
    countsB = zeros(nB,1,'logical');
else
    countsB = zeros(nB,1);
end

for i = 1:nB
    maskKey = ismember( A(:, keyVars), B(i, keyVars), 'rows' );
    if ~isempty(condFcn)
        Asub    = A(maskKey, :);
        maskKey = condFcn( Asub, B(i,:) );
    end
    countsB(i) = sum(maskKey);
end

end
