function fig_sample_species_composition()

ABUNDANCE_CUTOFF = 1; % minimum abundance to plot (percent)
SPECIES = ["Escherichia", "Citrobacter", "Enterobacter", "Pseudomonas", "Morganella", "Other", "Klebsiella"];
SPECIES_KEYS = lower(SPECIES);
n_species = numel(SPECIES);
other_idx = find(SPECIES == "Other");
species_idx = setdiff(1:n_species, other_idx);

add_required_paths

% Colors
cmixed = [0.8 0.8 0.8];
cmap = get_colors(SPECIES_KEYS);

% Load data
metaphlanAllT = loadMetaphlanAll();
metaphlanAllT = metaphlanAllT(metaphlanAllT.is_faecal & ~metaphlanAllT.is_rep, :);

% keep only genus-level rows (genus present, deeper ranks empty)
is_genus_level = metaphlanAllT.genus ~= "" & metaphlanAllT.species == "" & metaphlanAllT.sgb == "";
metaphlanGenusT = metaphlanAllT(is_genus_level, :);
metaphlanGenusT.genus = lower(metaphlanGenusT.genus);

toplotT = metaphlanGenusT;

% Get the most abundant genera
genus_abundance_summary = groupsummary(toplotT, 'genus', 'sum', 'abundance');
genus_abundance_summary = sortrows(genus_abundance_summary, 'sum_abundance', 'descend');
top_n = min(n_species - 1, height(genus_abundance_summary));
top_genera = genus_abundance_summary.genus(1:top_n);

if ~isempty(setdiff(top_genera, SPECIES_KEYS))
    warning('Top genera does not match species');
    top_genera
    SPECIES_KEYS
end

% create abundance matrix
fids = unique(toplotT.fid);
[~, j_fid] = ismember(toplotT.fid, fids);
[~, j_species] = ismember(toplotT.genus, SPECIES_KEYS);
j_species(j_species == 0) = other_idx;

record_value('Number of patients, metaphlan','num_pats_metaphlan',length(fids));

abundances = accumarray([j_fid, j_species], toplotT.abundance, [length(fids), n_species], @sum, 0);

assert(all(abs(sum(abundances, 2) - 100) <= 0.2));
abundances(:, other_idx) = 100 - sum(abundances(:, species_idx), 2);


[abundances, sids] = sortrows(abundances,[1 size(abundances,2)-1]);
fids = fids(sids);

is_mixed = ~any(abundances > 100-ABUNDANCE_CUTOFF, 2);
record_value('Number of mixed-genera patients, metaphlan','num_mixed_genera_pats',sum(is_mixed));

ecoli_col = strcmp(SPECIES_KEYS, 'escherichia');
n_ecoli = sum(abundances(:, ecoli_col) > 0);
record_value('Patients with Escherichia (%)','pct_pats_escherichia', n_ecoli / length(fids) * 100, '%.1f');

%% Figure 2b: genera abundance bar chart
figure(1); clf
ax = axes();
hbar = bar(ax, abundances(is_mixed,:), 'stacked', 'FaceColor','flat');
for k = 1:numel(hbar)
    hbar(k).FaceColor = cmap(k,:);
end

xlabel('Faecal pathometagenomes containing multiple genera')
ylabel('Genus abundance')
yticks = 0:20:100;
set(ax, ...
    'Position',[0.2 0.2 0.75 0.75],...
    'YLim',[0 100],...
    'fontsize',10,...
    'XTick',1:sum(is_mixed),...
    'XTickLabel', arrayfun(@(x) sprintf('P%i',x),fids(is_mixed), 'uniformoutput',false),...
    'YTick',yticks,...
    'YTickLabel',arrayfun(@(x) sprintf('%i%%',x),yticks,'UniformOutput',false))
legend(ax, SPECIES, 'location','se')

set(gca, 'fontsize',5, 'FontName','Helvetica')
set(gca, 'XTickLabelRotation',60)
set_axes_physical_size(gca, 4.1, 2, 'inch');
print_figure(gcf, 'fig_Genera_per_FITs')

%% Figure 2a: single-species vs multiple-species piechart
is_pure = ~is_mixed;
[~, dominant_genus_idx] = max(abundances(:, 1:numel(SPECIES)), [], 2);
pure_counts = accumarray(dominant_genus_idx(is_pure), 1, [numel(SPECIES), 1], @sum, 0);
present_pure_idx = find(pure_counts > 0);
mixed_count = sum(is_mixed);

toplot = pure_counts(present_pure_idx)';
labels = SPECIES(present_pure_idx);
slice_colors = cmap(present_pure_idx, :);

if mixed_count > 0
    toplot = [toplot, mixed_count];
    labels = [labels, "Multiple genera"];
    slice_colors = [slice_colors; cmixed];
end

figure(2); clf
frac = toplot / sum(toplot);
mixed_idx = find(labels == "Multiple genera", 1);
if isempty(mixed_idx)
    mixed_idx = numel(labels);
end
cum_before = sum(frac(1:mixed_idx-1));
start_angle = 360 * frac(mixed_idx) / 2;

p = pie(toplot);
for k = 1:numel(toplot)
    p(2*k-1).FaceColor = slice_colors(k,:);
    p(2*k).String = labels(k) + newline + "n=" + string(toplot(k));
    p(2*k).FontSize = 5;
    p(2*k).FontName = 'Helvetica';
end
view([90-start_angle, 90]);

% Print figure
set(gcf,'PaperPosition',[0 0 4 2])
set(gcf,'Color','none');
print_figure(gcf, 'fig_Genera_pie')

end
