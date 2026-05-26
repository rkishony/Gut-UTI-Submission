%% Community phenotyping colony count analysis
close all; clear; clc;
add_required_paths();

data_fldr = fullfile(get_folders('source_data'), 'Phenotyping');

exp_names = {'CommunityPhenotyping', 'CommunityPhenotypingB'};
MAX_COUNT = 50;

%% Load data and collect conditions
all_cond_labels = {};
all_counts = {};
all_mc_counts = {};
all_exp_labels = {};

for iExp = 1:numel(exp_names)
    fname = fullfile(data_fldr, exp_names{iExp}, 'results.csv');
    T = readtable(fname, 'VariableNamingRule', 'preserve');

    % Find NumColonies columns
    cols = T.Properties.VariableNames;
    nc_cols = cols(startsWith(cols, 'AllAbs_NumColonies_'));
    cond_names = cellfun(@(x) x(length('AllAbs_NumColonies_')+1:end), nc_cols, 'UniformOutput', false);

    % Filter out negative controls
    mask = ~T.is_neg_ctrl;

    note_cols = cols(startsWith(cols, 'AllAbs_Note_'));
    note_names = cellfun(@(x) x(length('AllAbs_Note_')+1:end), note_cols, 'UniformOutput', false);

    mc_idx = find(strcmp(cond_names, 'MC'));
    mc_counts = T{mask, nc_cols{mc_idx}};
    mc_note_idx = find(strcmp(note_names, 'MC'));
    mc_notes = T{mask, note_cols{mc_note_idx}};
    mc_counts(strcmpi(mc_notes, 'Lawn')) = MAX_COUNT;
    mc_counts = min(mc_counts, MAX_COUNT);

    for iC = 1:numel(cond_names)
        if strcmp(cond_names{iC}, 'MC')
            continue
        end
        short_exp = regexprep(exp_names{iExp}, 'CommunityPhenotyping', 'CP');
        if isempty(short_exp); short_exp = 'CP'; end
        all_cond_labels{end+1} = [short_exp ': ' cond_names{iC}]; %#ok
        counts = T{mask, nc_cols{iC}};
        ni = find(strcmp(note_names, cond_names{iC}));
        notes = T{mask, note_cols{ni}};
        counts(strcmpi(notes, 'Lawn')) = MAX_COUNT;
        counts = min(counts, MAX_COUNT);
        all_counts{end+1} = counts; %#ok
        all_mc_counts{end+1} = mc_counts; %#ok
        all_exp_labels{end+1} = exp_names{iExp}; %#ok
    end
end

nConds = numel(all_cond_labels);

%% Load per-FOBT_ID tables for cross-experiment comparison
for iExp = 1:numel(exp_names)
    fname = fullfile(data_fldr, exp_names{iExp}, 'results.csv');
    T = readtable(fname, 'VariableNamingRule', 'preserve');
    T = T(~T.is_neg_ctrl & ~isnan(T.FOBT_ID), :);

    cols = T.Properties.VariableNames;
    nc_cols = cols(startsWith(cols, 'AllAbs_NumColonies_'));
    note_cols = cols(startsWith(cols, 'AllAbs_Note_'));
    cond_names = cellfun(@(x) x(length('AllAbs_NumColonies_')+1:end), nc_cols, 'UniformOutput', false);
    note_names = cellfun(@(x) x(length('AllAbs_Note_')+1:end), note_cols, 'UniformOutput', false);

    for iC = 1:numel(cond_names)
        cname = cond_names{iC};
        c = T{:, nc_cols{iC}};
        ni = find(strcmp(note_names, cname));
        n = T{:, note_cols{ni}};
        c(strcmpi(n, 'Lawn')) = MAX_COUNT;
        c = min(c, MAX_COUNT);
        exp_data(iExp).counts.(strrep(cname,'-','_')) = c; %#ok
    end
    exp_data(iExp).FOBT_ID = T.FOBT_ID; %#ok
end

%% Cross-experiment agreement figure
THRESH = 25;
pairs = {
    'FOS-16', 'FOS', 'CPB:FOS-16 vs CP:FOS'
    'FOS-32', 'FOS', 'CPB:FOS-32 vs CP:FOS'
    'NIT-8',  'NIT', 'CPB:NIT-8 vs CP:NIT'
};
nPairs = size(pairs, 1);

[~, ia, ib] = intersect(exp_data(2).FOBT_ID, exp_data(1).FOBT_ID);
n_cp_only  = sum(~ismember(exp_data(1).FOBT_ID, exp_data(2).FOBT_ID));
n_cpb_only = sum(~ismember(exp_data(2).FOBT_ID, exp_data(1).FOBT_ID));
fprintf('Matched samples: %d\n', numel(ib));
fprintf('CP only (not in CPB): %d\n', n_cp_only);
fprintf('CPB only (not in CP): %d\n', n_cpb_only);


fig2 = figure('Position', [100 100 600 250*nPairs]);
for iP = 1:nPairs
    cpb_field = strrep(pairs{iP,1}, '-', '_');
    cp_field  = strrep(pairs{iP,2}, '-', '_');
    c_cpb = exp_data(2).counts.(cpb_field)(ia);
    c_cp  = exp_data(1).counts.(cp_field)(ib);

    % Scatter
    subplot(nPairs, 2, (iP-1)*2 + 1);
    scatter(c_cp, c_cpb, 10, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    plot([0 MAX_COUNT], [0 MAX_COUNT], 'k--', 'LineWidth', 0.5);
    xlim([0 MAX_COUNT]); ylim([0 MAX_COUNT]);
    xlabel(['CP: ' pairs{iP,2}]); ylabel(['CPB: ' pairs{iP,1}]);
    title(pairs{iP,3}, 'Interpreter', 'none');

    % Confusion matrix (threshold = THRESH)
    subplot(nPairs, 2, (iP-1)*2 + 2);
    r_cp  = c_cp  >= THRESH;
    r_cpb = c_cpb >= THRESH;
    cm = [sum(~r_cp & ~r_cpb), sum(r_cp & ~r_cpb);
          sum(~r_cp &  r_cpb), sum(r_cp &  r_cpb)];
    imagesc(cm); colormap(gca, flipud(gray));
    set(gca, 'XTick', [1 2], 'XTickLabel', {'<25','>=25'}, ...
             'YTick', [1 2], 'YTickLabel', {'<25','>=25'});
    xlabel(['CP: ' pairs{iP,2}]); ylabel(['CPB: ' pairs{iP,1}]);
    title(sprintf('Agreement: %.0f%%', 100*sum(r_cp==r_cpb)/numel(r_cp)));
    for r = 1:2
        for c = 1:2
            text(c, r, num2str(cm(r,c)), 'HorizontalAlignment','center', ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end

%% Plot (Figure 1)
count_edges = 0:2:MAX_COUNT;
norm_edges = 0:0.05:2;

fig = figure('Position', [100 100 900 120*nConds + 80]);
axs = gobjects(nConds, 3);

for iC = 1:nConds
    counts = all_counts{iC};
    mc = all_mc_counts{iC};
    norm_counts = counts ./ max(mc, 1);
    norm_counts = min(norm_counts, norm_edges(end));
    is_bottom = (iC == nConds);

    % Column 1: histogram of counts
    axs(iC,1) = subplot(nConds, 3, (iC-1)*3 + 1);
    histogram(counts, count_edges);
    xlim([0 MAX_COUNT]);
    ylabel(all_cond_labels{iC}, 'Interpreter', 'none');
    if is_bottom
        xlabel('Colony count');
        xticks([0 10 20 30 40 50]);
        xticklabels({'0','10','20','30','40','50+'});
    else
        set(gca, 'XTickLabel', []);
    end

    % Column 2: histogram of normalized counts
    axs(iC,2) = subplot(nConds, 3, (iC-1)*3 + 2);
    histogram(norm_counts, norm_edges);
    xlim([0 norm_edges(end)]);
    if is_bottom
        xlabel('Count / MC count');
    else
        set(gca, 'XTickLabel', []);
    end

    % Column 3: scatter count vs MC count
    axs(iC,3) = subplot(nConds, 3, (iC-1)*3 + 3);
    scatter(mc, counts, 10, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    plot([0 MAX_COUNT], [0 MAX_COUNT], 'k--', 'LineWidth', 0.5);
    xlim([0 MAX_COUNT]); ylim([0 MAX_COUNT]);
    if is_bottom
        xlabel('MC count');
        xticks([0 10 20 30 40 50]);
        xticklabels({'0','10','20','30','40','50+'});
    else
        set(gca, 'XTickLabel', []);
    end
    yticks([0 10 20 30 40 50]);
    yticklabels({'0','10','20','30','40','50+'});
end

% Squeeze vertical gaps
set(axs, 'FontSize', 8);
for iC = 1:nConds
    for jC = 1:3
        pos = get(axs(iC,jC), 'Position');
        vgap = 0.02;
        h = (0.88 - vgap*(nConds-1)) / nConds;
        top = 0.95 - (iC-1)*(h + vgap);
        set(axs(iC,jC), 'Position', [pos(1), top - h, pos(3), h]);
    end
end
