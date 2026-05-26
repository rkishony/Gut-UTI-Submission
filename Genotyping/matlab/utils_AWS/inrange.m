function tf = inrange(x, ab, includeEndpoints)
% INRANGE   True where a < x < b (or a a< x a< b if includeEndpoints=true)
if nargin<4, includeEndpoints = true; end
a = ab(1);
b = ab(2);
if includeEndpoints
    tf = (x >= a) & (x <= b);
else
    tf = (x >  a) & (x <  b);
end
end
