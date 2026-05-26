function A = load_alignment(matched, include_rep)
if nargin<1
    matched = true;
end
if nargin<2
    include_rep = false;
end

folders = get_folders();
pthA = fullfile(folders.from_labserver, 'summary', 'alignment_stats_coalesced.csv');

A = readtable(pthA);
T = load_trimstats();

[~, j] = ismember(A.uti_run_name, T.run_name);
A.uti_sample_name = T.sample_name(j);
A.uti_is_rep = T.is_rep(j);
A.UTI_ID = T.uid(j);
A.UTI_FID = T.fid(j);
A.UTI_isolate = T.is_isolate(j);

A.etha0 = A.uncov_len ./ A.genome_len;

if matched
    A = A(strcmp(A.pair_type, 'matched'), :);
end

if ~include_rep
    A = A(~A.uti_is_rep, :);
end

end

