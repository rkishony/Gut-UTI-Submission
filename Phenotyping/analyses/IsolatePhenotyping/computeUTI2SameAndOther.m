function [D2Same, vargout] = computeUTI2SameAndOther(j_resistance_maps, j_fobt, j_uti, j_isfobt, j_isuti, uti_ids, uti_patient_ids, compute_d2other_too)
patients = unique(uti_patient_ids);
count_utis = 0;
D2Same = [];
D2Other = [];

j_uti(isnan(j_uti))=0;
j_isuti(isnan(j_isuti))=0;

for u = 1:length(uti_ids)
    this_u = uti_ids(u);
    uI = find(j_uti==this_u & j_isuti);
    
    this_f = uti_patient_ids(u);

    assert(length(this_f)==1)
    
    if this_f==0
        fprintf('no FOBT isolate for UTI %d.\n')
        continue % is never reached, sanity check
    end
    if isempty(find(j_uti==this_u & j_isuti)) 
        fprintf('no UTI isolate grew without AB for UTI %d.\n',this_u);
        continue % is never reached, sanity check
    end
    if isempty(find(j_fobt==this_f & j_isfobt))
        continue % this condition was added to exclude microbiomes without a UTI in the dataset but is not necessarily important to keep
    end
    
    fI = find(j_fobt==this_f & j_isfobt);
    
    sUFdist = computeUFdist(j_resistance_maps, uI, fI);
    D2Same = [D2Same; sUFdist];
    
    clear fI
    clear sUFdist
    if compute_d2other_too   
        these_other = [];
        for f = 1:length(patients)
            other_f = patients(f);
            if other_f~=this_f
                fI = find(j_fobt==other_f & j_isfobt);
                
                oUFdist = computeUFdist(j_resistance_maps, uI, fI);
                
                these_other = [these_other; oUFdist];
                clear oUFdist
            end
        end
        D2Other = [D2Other; these_other];
    end
    clear uI
    count_utis = count_utis+1;
end
vargout = D2Other;

    function UFdist = computeUFdist(j_resistance_maps, uI, fI)
        this_u_profs = j_resistance_maps(uI,:);
        this_f_profs = j_resistance_maps(fI,:);
        UFdists = pdist2(this_u_profs,this_f_profs,@distfun);
        UFdist = min(UFdists, [], 'all');
    end

    function d2 = distfun(ZI, ZJ)
        notnan = ~isnan(ZI) & ~isnan(ZJ);   % positions where both are numbers
        diffs = abs(ZI - ZJ);               % compute differences
        diffs(~notnan) = 0;                 % ignore NaNs by zeroing them out
        d2 = sum(diffs, 2);                 % sum across columns
    end
end