function fig_MLST_heatmap_and_bar()
% figure 2c heatmap of isolate MLST

% species identifiers (as they appear in the MLST table)
species = ["ecoli", "klebsiella", "cfreundii", "ecloacae"];
panel_species = [species, "other"];

species_nice_names = containers.Map(...
    {'ecoli',   'klebsiella',    'cfreundii',   'ecloacae',    'other'}, ...
    {'E. coli', 'K. pneumoniae', 'C. freundii', 'E. cloacae',  'Other'});

% map species identifiers to genus-level color keys in colors.txt
species_to_genus = containers.Map(...
    {'ecoli',       'klebsiella', 'cfreundii',  'ecloacae',     'other'}, ...
    {'escherichia', 'klebsiella', 'citrobacter', 'enterobacter', 'other'});

RARE_ST_THRESHOLD = 2;  % STs less than this many occurrences are collapsed into "other" bucket
MIN_DISPLAY_FRACTION = 0.15;  % floor for non-zero cells so all present STs are visually distinguishable

% load data
add_required_paths();
mlstT = loadisolateMLST();
mlstT = mlstT(mlstT.is_faecal, :);
is_ok = mlstT.status == "ok";

% stats
stat = struct;
stat.tot = height(mlstT);
stat.with_mlst = sum(is_ok);
stat.no_assembly = sum(mlstT.status == "missing_file");
stat.no_call = sum(mlstT.status == "no_call");
stat.num_pats = length(unique(mlstT.fid(is_ok)));
record_value('Number of fecal isolates', 'num_faecal_iso', stat);

mlstT = mlstT(is_ok, :);

% Unknown STs -> "other" (0) for display/collapsing
assert(sum(mlstT.ST == 0) == 0);
is_unknown_ST = isnan(mlstT.ST);
mlstT.ST(is_unknown_ST) = 0;
mlstT.ST_label = arrayfun(@(x) sprintf('ST%d', x), mlstT.ST, 'UniformOutput', false);
mlstT.ST_label(is_unknown_ST) = {'unknown'};

% Count occurrences of each (fid, mlst) pair
[~, i_mlsts, j_mlsts] = unique(mlstT(:, {'species', 'ST'}), 'rows');
mlsts = mlstT(i_mlsts, :);
[patients, ~, j_patients] = unique(mlstT(:, {'fid'}));
counts = accumarray([j_patients, j_mlsts], 1, [height(patients), height(mlsts)]);
mlsts.is_rare = (sum(counts, 1)' < RARE_ST_THRESHOLD);
mlsts.is_other = mlsts.ST == 0;

% Sort patients and mlsts by number of unique MLSTs per patient
patients.mlst_num = sum(counts > 0, 2);
mlsts.patient_num = sum(counts > 0, 1)';
[patients, spat] = sortrows(patients, 'mlst_num', 'descend');
[mlsts, smlst] = sortrows(mlsts, 'patient_num', 'descend');
counts = counts(spat, smlst);

% Normalize by number of isolates per patient
counts = counts ./ sum(counts, 2);

figure(7)
set(gcf, 'Visible', 'off')
clf

gap_h = 0.01;
cbwidth = 0.01;
cbheight = 0.4;
ylabelswidth = 0.05;
margin_h = 0.12;
margin_b = 0.2;
margin_t = 0.1;
hmfraction = 0.5;
hmheight = 1 - margin_t - margin_b;
bpfraction = 1 - hmfraction - cbwidth - ylabelswidth - 2*margin_h;

% colors
genus_keys = cellfun(@(s) species_to_genus(s), cellstr(panel_species), 'UniformOutput', false);
species_colors = get_colors(genus_keys);
n_cmap = 128;
nonzero_shade = linspace(0.1, 1, n_cmap)' .^ 0.7;
nonzero_gray = 1 - nonzero_shade;
cmap = [1 1 1; repmat(nonzero_gray, 1, 3)];

num_species = length(species);
num_panels = length(panel_species);
num_cols_per_species = zeros(num_panels, 1);
for i = 1:num_species
    j_mlsts = find(mlsts.species == species(i));
    is_other_bucket = mlsts.is_rare(j_mlsts) | mlsts.is_other(j_mlsts);
    num_cols_per_species(i) = sum(~is_other_bucket) + any(is_other_bucket);
end
num_cols_per_species(num_species+1) = 1;  % "other" species panel with single column
total_plot_cols = sum(num_cols_per_species);
widths = (num_cols_per_species / total_plot_cols) * (hmfraction - num_panels*gap_h);

first_panel_axis = [];
for i = 1:num_panels
    if i <= num_species
        j_mlsts = find(mlsts.species == species(i));
        mlsts_i = mlsts(j_mlsts, :);
        counts_i = counts(:, j_mlsts);

        is_other = mlsts_i.is_rare | mlsts_i.is_other;
        plot_counts_i = counts_i(:, ~is_other);
        labels_i = mlsts_i.ST_label(~is_other)';
        if any(is_other)
            plot_counts_i = [plot_counts_i, sum(counts_i(:, is_other), 2)];
            labels_i = [labels_i, {'others'}];
            n_others_per_pat = sum(counts_i(:, is_other) > 0, 2);
        else
            n_others_per_pat = [];
        end
    else
        is_other_species = ~ismember(mlsts.species, species);
        plot_counts_i = sum(counts(:, is_other_species), 2);
        labels_i = {'any'};
        n_others_per_pat = sum(counts(:, is_other_species) > 0, 2);
    end
    n_cols = size(plot_counts_i, 2);
    assert(n_cols == num_cols_per_species(i));
    left = margin_h + ylabelswidth + cbwidth + sum(widths(1:i-1)) + (i-1)*gap_h;
    h1 = axes('Position',[left, margin_b, widths(i), hmheight]);
    plot_display_i = plot_counts_i;
    % plot_display_i(plot_display_i > 0) = max(plot_display_i(plot_display_i > 0), MIN_DISPLAY_FRACTION);
    imagesc(h1, plot_display_i)
    set(h1,'fontsize',5,...
        'ytick',1:height(patients), 'yticklabel',[],...
        'box','on',...
        'ydir','reverse',...
        'ylim',[0.4 height(patients)+0.6],...
        'xlim',[0.4 n_cols+0.6],'xtick',1:n_cols,...
        'XTickLabel',labels_i,'XTickLabelRotation',90,...
        'colormap',cmap);
    h1.XAxis.TickLabelGapOffset = 2;
    if i==1
        set(h1, 'YTickLabel', arrayfun(@(x)sprintf('P%i',x), patients.fid, 'UniformOutput', false));
        h1.YAxis.TickLabelGapOffset = 2;
        ylabel('Patients','FontSize',6);
        first_panel_axis = h1;
    end
    set(h1, 'CLim', [0 1]);
    % Annotate "others"/"any" cells where >1 ST is collapsed into the bucket
    % for r = find(n_others_per_pat > 1)'
    %     text(h1, n_cols, r, num2str(n_others_per_pat(r)), ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    %         'FontSize', 4, 'Color', 'k');
    % end
    xlabel(h1, species_nice_names(char(panel_species(i))), 'FontSize', 5, 'FontAngle', 'italic');

    xlim_ = xlim(gca); ylim_ = ylim(gca);
    rectangle(gca, 'Position', [xlim_(1), ylim_(1), diff(xlim_), diff(ylim_)], ...
              'EdgeColor', species_colors(i, :), ...
              'LineWidth', 2, ...
              'LineStyle', '-', 'Clipping','off')
end

% colorbar
c1 = colorbar(first_panel_axis, 'Location', 'westoutside', 'ytick', 0:0.1:1);
c1.Position = [margin_h*0.8 margin_b+cbwidth-gap_h+hmheight-cbheight cbwidth cbheight];
c1.Label.String = {'Fraction of isolates', 'from faecal pathobiome'};
c1.Label.FontSize = 6;
c1.FontSize = 5;

% bar plot of number of unique MLSTs per patient
h2 = axes('position',[1-bpfraction-margin_h, margin_b, bpfraction, 1-margin_t-margin_b]);
mlst_num = patients.mlst_num;

barh(1:size(counts,1), mlst_num, 'barwidth', 1, 'FaceColor', [0.5 0.5 0.5])
set(gca,'YTickLabel',[],'YDir','reverse','FontSize',6,...
    'ylim',[0.5 size(counts,1)+0.5],'XTick',0:max(mlst_num))
xlim(gca, [0 max(mlst_num)+0.5]);

xlabel('Number of unique MLSTs per patient','FontSize',6)
record_value('Number of patients with >1 ST','num_pats_w_multi_MLST',sum(mlst_num>1));
record_value('Average number of unique STs per patient','avg_MLSTs_per_pat',mean(mlst_num),'%.2f');

% export figure
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 7.5 2.25]);
print_figure(gcf, 'fig_MLST_heatmap_and_bar')

end