function T = load_trimstats

folders = get_folders();
pthT = fullfile(folders.from_labserver, 'summary', 'trim_stats.csv');

T = readtable(pthT);
T.is_rep = strcmp(T.is_rep, 'True');

T.is_isolate = strcmp(T.is_isolate, 'True');
T.is_faecal = strcmp(T.is_faecal, 'True');
T.has_two_fasta = strcmp(T.has_two_fasta, 'True');
T.downsampled = strcmp(T.downsampled, 'True');

end

