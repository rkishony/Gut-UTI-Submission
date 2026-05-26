function etha0s_and_distances_data = analyze_relation_between_UTIcovbre_and_FITprofiles(etha0_uti_gid_fobt_gid, MHSdata, personalmats, patients)

% minimal distance analysis
all_info = getMinDist_MHS(etha0_uti_gid_fobt_gid, MHSdata, personalmats, patients);

% cap etha0s
capped_etha0s = all_info.etha0;
capped_etha0s(round(capped_etha0s,4) <= 1e-4) = 1e-4;
f_capped_a_etha0s = log10(capped_etha0s); % used to be log10(1-breadths)

% put capped etha0 into table with everything else and get number of
% isolates per patients into output table
etha0s_and_distances_data = all_info;
etha0s_and_distances_data.capped_etha0 = f_capped_a_etha0s;
n_isolates = arrayfun(@(c) size(all_info.fobt_isolate_profiles{c},1), 1:height(all_info));
etha0s_and_distances_data.n_isolates = n_isolates';


end