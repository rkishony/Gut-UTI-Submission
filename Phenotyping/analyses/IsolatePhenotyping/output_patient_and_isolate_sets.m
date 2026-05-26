clear

all_res = readtable(fullfile(get_folders('source_data'), 'Phenotyping', 'IsolateResData.csv'), 'TextType', 'string');

fprintf('There were a total of %d sets of 96 wells plates (though some were arrayed as 384).\n', height(all_res)/96)


%% get patient and isolate numbers from every experiment/used in every analysis:

S = summarize_isolates_in_experiments(all_res);

A = summarize_isolates_in_analyses(all_res);
