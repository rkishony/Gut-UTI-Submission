function [d, nMatch, matchRatio] = nansDistance(x, y, minOverlap)
    if nargin < 3
        minOverlap = 1;
    end

    valid = ~isnan(x) & ~isnan(y);
    nMatch = sum(valid);

    % number of positions that could possibly match (at least one is non-NaN)
    possible = sum(~isnan(x) | ~isnan(y));

    if possible == 0
        matchRatio = 0;
    else
        matchRatio = nMatch / possible;
    end

    if nMatch < minOverlap
        d = inf;
        return
    end

    % normalized distance (recommended)
    d = norm(x(valid) - y(valid),1);% / sqrt(nMatch);
end