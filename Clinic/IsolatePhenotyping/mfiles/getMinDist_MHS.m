function results = getMinDist_MHS(etha0_vars, MHSdata, pooleddata, patients, varargin)

% 'patients' is the labels of 'pooleddata'

etha_0 = etha0_vars(:,1);  uti_id = etha0_vars(:,2);  fobt_id = etha0_vars(:,3);

rel_patients = intersect(patients, fobt_id);

% must have an experimentally measured profile to be in comparison, taking
% only those from MHSdata:
MHSdata = MHSdata(ismember(MHSdata.FOBT_ID,rel_patients),:);

up = rel_patients; 

results = table('Size',[0 9],'VariableTypes',{'double','double','cell','cell','double','cell','double','double','double'}, ...
    'VariableNames',{'UTI_ID','patientID','UTI_profile','UTI_prof2compare','etha0','fobt_isolate_profiles','distance','nMatches','MatchRatio'});

count = 0;
% loop over patients:
for pa = up'
    
    if ismember(pa,patients)
        
        ug = uti_id(fobt_id==pa);
        bg = etha_0(fobt_id==pa);

        cur_fobt_profs = pooleddata{patients==pa};
        if isempty(cur_fobt_profs)
            sprintf('no fobt isolates for patient %d', pa)
            continue
        end

        % loop over this patient's UTI IDs
        for i = 1:numel(ug)
            
            count = count+1;

            current_UTI = ug(i);
            UTI_ID = current_UTI;
            
            % find MHS profiles of this UTI ID
            iii = find(MHSdata.UTI_ID==current_UTI);
            if isempty(iii)
                continue
            else
                current_profile = MHSdata.vnum(iii,:);
            end
            
            patientID = pa;
            UTI_profile = {current_profile};
            
            % keep nans but take only 'R' to be resistant
            current_profile(current_profile <2)  = 0;
            current_profile(current_profile==2) = 1;
            UTI_prof2compare = {current_profile};

            current_etha0 = bg(i);
            e0 = current_etha0;
            fobt_isolate_profiles = {cur_fobt_profs};

            % if more than one MHS profile is available with Maccabi,
            % choose the one with more resistance measurements:

            notmeas = isnan(current_profile);

            [measnum, more_measured_one] = min(sum(notmeas,2));

            % if only one UTI has the most measurements, take this one.
            % otherwise loop over them one by one and take the pair
            % with the minimal distance.
            if sum(sum(notmeas,2)==measnum)>1
                current_profile = current_profile(more_measured_one,:);
                notmeas = isnan(current_profile);
            end
            
            % find minimal distance: for each UTI culture, compare profile(s) to all FOBT
            % measurements
            bestDist = inf;
            bestNMatch = nan;
            bestMatchRatio = nan;
            for u = 1:size(current_profile,1)
                cur_uprof = current_profile(u,:);
                for fi = 1:size(cur_fobt_profs,1)
                    cur_fprof = cur_fobt_profs(fi,:);

                    [d, nMatch, matchRatio] = nansDistance(cur_uprof, cur_fprof, 2);
                    if d<bestDist
                        bestDist = d;
                        bestNMatch = nMatch;
                        bestMatchRatio = matchRatio;
                    end
                end
            end
            distance = bestDist;
            nMatches = bestNMatch;
            MatchRatio = bestMatchRatio;
            newrow = {UTI_ID, patientID, UTI_profile, UTI_prof2compare, e0, fobt_isolate_profiles, distance, nMatches, MatchRatio};
            results(end+1,:) = newrow;

        end
    end
end

end
