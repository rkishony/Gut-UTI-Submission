function create_all_figures()

add_required_paths();
record_value([])

figure_scripts = {
    'fig_feacal_uti_cov'
    'fig_MLST_heatmap_and_bar'
    'fig_sample_species_composition'
};

for i = 1:numel(figure_scripts)
    fprintf('\n%s\nRunning %s\n%s\n', repmat('=',1,60), figure_scripts{i}, repmat('=',1,60));
    try
        feval(figure_scripts{i});
    catch e
        fprintf('WARNING: %s failed: %s\n', figure_scripts{i}, e.message);
    end
    close all;
end

records = record_value();
writetable(records, fullfile(get_folders('figures'), 'record_values.csv'));
end
