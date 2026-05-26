function mV4 = convert_v123count_to_v4(C123, count_type_dim, time_method, time_dim)
% Conver an array of resistance counts to array of resistance status.
%
% C123           : Input array of UTI counts where dimension `dim_res` of size 3,
%                  represent count-type: 
%                  (1) Null measurements, (2) Sensitive, (3) Resistant
%
% count_type_dim : The count-type dimension. Last dimension if not provided.
%
% time_method    : Indicates how PUC measurements are aggregated over time:
%                  'local' : Each bin is treated separately (default).
%                  'all'   : The results of the bin and all future bins.
%                  'first' : The first UTI from this bin forward.
%                  'first_measured'
%                          : The first measurement from this bin forward.
%
% time_dim       : The array of timebins. 
%                  Required only if method is not 'local'.
%
% OUTPUT
% mV4  : an array of same size as C123 except singleton at dim.
%        Colapsing this dimension, the output vector values are:
%
%        1 - NONE No PUCs.
%        2 - NULL Only PUCs with no measurements.
%        3 - SEN  At least 1 sensitive PUC and no resistant PUCs.
%        4 - RES  At least 1 resistant PUC.
%

if nargin<2 || isempty(count_type_dim)
    count_type_dim = ndims(C123);
end
if nargin<3
    time_method = 'local';
end

n = 3;
assert(size(C123, count_type_dim) == n)

% Convert counts to 1-4 value:
sz = size(C123);
otherdims = setdiff(1:numel(sz),count_type_dim);
perm = [count_type_dim, otherdims];
m = numel(C123) / n;
C123 = permute(C123, perm);

C123 = reshape(C123, [n, m]);

is_puc      = sum(C123, 1) > 0;
is_measured = sum(C123(2:3, :), 1) > 0;
is_res      = C123(3,:) > 0;

V4 = 1 + is_puc + is_measured + is_res;
V4 = reshape(V4, [1, sz(otherdims)]);
V4 = ipermute(V4, perm);

% Aggregate timebins:

switch time_method
    case 'local'
        mV4 = V4;
    case 'all'
        mV4 = apply_backwards(@cummax, V4, time_dim);
    case {'first', 'first_measured'}
        V4_first = get_next_marked_value(V4, V4>1, time_dim);
        V4_first(isnan(V4_first)) = 1;
        if strcmp(time_method,'first')
            mV4 = V4_first;
        else
            V4_first_measured = get_next_marked_value(V4, V4>2, time_dim);
            V4_first_measured(isnan(V4_first_measured)) = V4_first(isnan(V4_first_measured));
            mV4 = V4_first_measured;
        end
    otherwise
        error('unknown time_method %s', time_method)
end

end
