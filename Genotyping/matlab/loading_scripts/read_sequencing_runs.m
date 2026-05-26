function T = read_sequencing_runs(csvPath)

if nargin < 1
    csvPath = fullfile(get_folders('analyses_server'), 'summary', 'trim_stats.csv');
end

T = readtable(csvPath, 'TextType', 'string');

boolCols = ["is_isolate", "is_faecal", "is_rep", "has_two_fasta", "downsampled"];
for c = boolCols
    if ismember(c, T.Properties.VariableNames)
        v = T.(c);
        is_true = strcmpi(v, "True");
        is_false = strcmpi(v, "False");
        assert(all(xor(is_true, is_false)), 'Column %s contains non-boolean values', c);
        T.(c) = is_true;
    end
end

end
