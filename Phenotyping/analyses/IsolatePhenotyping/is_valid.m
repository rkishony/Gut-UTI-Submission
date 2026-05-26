function isValid = is_valid(T, isMeasured, n)
if nargin<3
    n=0; % nan = no growth in MC. 
end

res = get_resistance_matrix(T);

isValid = isMeasured & ~(sum(isnan(res),2)>n); 
end