function [hbin, binedges] = breadths_and_distances_into_colorcoded_bars_data(breadths_and_distances_data, r)

etha0 = breadths_and_distances_data.capped_etha0;
dist = breadths_and_distances_data.distance;
UTIid = breadths_and_distances_data.UTI_ID;
nisolates = breadths_and_distances_data.n_isolates;

binedges = -4.25:0.5:0.25;
for e = 1:length(binedges)-1
    i_hbins{e} = (binedges(e) < etha0) & (binedges(e+1) >= etha0);
end

withinrange = nisolates >= r(1) & nisolates <= r(2); % if we want to filter by isolate number so that it doesn't vary so much, possibly affecting the results.
unchecked = UTIid;
for i = 1:length(i_hbins)
    identical = sum(i_hbins{i} & dist==0 & withinrange);
    diff1  = sum(i_hbins{i} & dist>0 & dist<=1 & withinrange);
    over1 = sum(i_hbins{i} & dist>1 & withinrange);
  
    hbin(i,:) = [identical, diff1, over1];

    utis_identical = UTIid(i_hbins{i} & dist==0 & withinrange);
    utis_diff1 = UTIid(i_hbins{i} & dist>0 & dist<=1 & withinrange);
    utis_over1 = UTIid(i_hbins{i} & dist>1 & withinrange);
    ubin{i,:} = {utis_identical, utis_diff1, utis_over1};

    %%%
    if ismember(utis_identical,UTIid) 
        unchecked(ismember(unchecked, utis_identical)) = 0;
    end
    if ismember(utis_diff1,UTIid)
        unchecked(ismember(unchecked, utis_diff1)) = 0;
    end
    if ismember(utis_over1,UTIid)
        unchecked(ismember(unchecked, utis_over1)) = 0;
    end
end
end