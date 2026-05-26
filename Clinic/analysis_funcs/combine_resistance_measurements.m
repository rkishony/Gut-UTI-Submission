function vnum = combine_resistance_measurements(vnum, method, dim)
% vnum: MHS array of nan, 0, 1, 2, 3 of multiple antibiotics across the
% `dim` dimension.
%
% return an array of one effective antibiotic combining all the raw
% resistances.

switch method
    case 'max'
        % resistant if resistant to least one drug
        vnum = max(vnum, [], dim, 'omitnan'); 
    case 'min'
        % sensitive if sensitive to least one drug
        vnum = min(vnum, [], dim, 'omitnan'); 
    case 'median'
        % median among valid measurements
        vnum = median(vnum, dim, 'omitnan');
    otherwise
        error('unrecognized method `%`', method)
end

end