function metaphlanAllT = loadMetaphlanAll()

strainphlan_path = fullfile(get_folders('analyses_server'), 'strainphlan');
profile_files = dir(fullfile(strainphlan_path, 'F*', 'profile.tsv'));

rows = cell(0, 10);
for ifile = 1:length(profile_files)
    profile_file = fullfile(profile_files(ifile).folder, profile_files(ifile).name);
    [~, run_name] = fileparts(profile_files(ifile).folder);

    fid = fopen(profile_file);

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        if isempty(line) || startsWith(line, '#')
            continue
        end

        line_columns = strsplit(line, sprintf('\t'));
        abundance = str2double(line_columns{3});

        rank_values = parseCladeName(line_columns{1});
        rank_cells = num2cell(rank_values);
        rows(end+1, :) = {run_name, rank_cells{:}, abundance};
    end
    fclose(fid);
end

metaphlanAllT = cell2table(rows, ...
    'VariableNames', {'run_name','kingdom','phylum','class','order','family','genus','species','sgb','abundance'});

sequencing_runsT = read_sequencing_runs();
[~, run_idx] = ismember(metaphlanAllT.run_name, sequencing_runsT.run_name);
assert(all(run_idx > 0), 'Some strainphlan runs have no sequencing-runs metadata');

metaphlanAllT.fid = sequencing_runsT.fid(run_idx);
metaphlanAllT.is_faecal = sequencing_runsT.is_faecal(run_idx);
metaphlanAllT.is_rep = sequencing_runsT.is_rep(run_idx);
end


function rank_values = parseCladeName(clade_name)
    tokens = strsplit(clade_name, '|');
    rank_values = strings(1, 8);
    ntok = min(length(tokens), 8);
    expected_rank_codes = 'kpcofgst';  % kingdom, phylum, class, order, family, genus, species, sgb
    for itok = 1:ntok
        token = tokens{itok};
        assert(length(token) >= 4);
        assert(strcmp(token(2:3), '__'));
        assert(token(1) == expected_rank_codes(itok));
        rank_values(itok) = string(token(4:end));
    end
end
