function [cov_same, cov_other] = fig_feacal_uti_cov()

%% Parameters
FOCAL_UTI_ISO   = 'U15A_hs1';
FOCAL_OTHER_FID = 109;
ETA0_CLAMP      = 1e-4;
AGGREGATION_LEVEL = 'isolate';  % 'isolate' | 'uti' | 'patient'
FIRST_ISOLATE_ONLY = true;      % true: keep only the first isolate (lowest letter) per UTI

csamepat = get_colors('FAECAL_same');
cdiffpat = get_colors('FAECAL_different');

%% Setup
add_required_paths();
output_dir = get_folders('analyses_server');

coalT = read_stats_table(output_dir, 'alignment_stats_coalesced.csv');
runT  = read_stats_table(output_dir, 'alignment_stats.csv');
trimT = read_sequencing_runs();

% Add uid, fid columns to coalT
focal_uti = trimT(trimT.run_name == FOCAL_UTI_ISO, :);
[~, uti_loc] = ismember(coalT.uti_run_name, trimT.run_name);
coalT.uti_uid = trimT.uid(uti_loc);
coalT.uti_fid = trimT.fid(uti_loc);
coalT.uti_isolate = trimT.isolate(uti_loc);
first_faecal_run = extractBefore(coalT.faecal_run_names + ",", ",");
[~, loc] = ismember(first_faecal_run, trimT.run_name);
coalT.faecal_fid = trimT.fid(loc);

if FIRST_ISOLATE_ONLY
    coalT = sortrows(coalT, 'uti_isolate');
    [~, ia] = unique(coalT(:, {'uti_uid', 'faecal_fid'}), 'rows', 'first');
    coalT = coalT(ia, :);
end

% Add faecal_fid to runT
[~, loc] = ismember(runT.faecal_run_name, trimT.run_name);
runT.faecal_fid = trimT.fid(loc);

%% ========================================================================
%  Figure 1: CDF of eta_0  (same-patient vs different-patient)
%  ========================================================================

ok = coalT.status ~= "no_log" & coalT.genome_len > 0;
coalT = coalT(ok, :);
coalT.eta0 = coalT.uncov_len ./ coalT.genome_len;

% Aggregate eta0 according to AGGREGATION_LEVEL
switch AGGREGATION_LEVEL
    case 'isolate'
        [~, ia, G] = unique(coalT(:, {'uti_run_name', 'faecal_fid'}), 'rows');
    case 'uti'
        [~, ia, G] = unique(coalT(:, {'uti_uid', 'faecal_fid'}), 'rows');
    case 'patient'
        [~, ia, G] = unique(coalT(:, {'uti_fid', 'faecal_fid'}), 'rows');
end
eta0_vec = accumarray(G, coalT.eta0, [], @min);

% focal points: find via G before clamping
focal_fids = [focal_uti.fid, FOCAL_OTHER_FID];
j_focal = zeros(size(focal_fids));
for k = 1:numel(focal_fids)
    jj = find(coalT.uti_run_name == FOCAL_UTI_ISO & coalT.faecal_fid == focal_fids(k));
    assert(isscalar(jj));
    if coalT.eta0(jj) ~= eta0_vec(G(jj))
        error('Focal isolate %s is not the @min of its F%d group', FOCAL_UTI_ISO, focal_fids(k));
    end
    j_focal(k) = jj;
end
eta0_focal_same = eta0_vec(G(j_focal(1)));
eta0_focal_other = eta0_vec(G(j_focal(2)));

% separate cdf for same- and different-patient
is_same_vec = coalT.pair_type(ia) == "matched";
eta0_same  = sort(eta0_vec(is_same_vec),  'ascend');
eta0_other = sort(eta0_vec(~is_same_vec), 'ascend');
p_same  = linspace(0, 1, numel(eta0_same));
p_other = linspace(0, 1, numel(eta0_other));

% Record stats
record_value('Number of same-patient sample pairs', 'num_same_pat_sample_pairs', numel(eta0_same));
record_value('Number of different-patient sample pairs', 'num_diff_pat_sample_pairs', numel(eta0_other));

[~, p_ks] = kstest2(eta0_same, eta0_other);
record_value('KS-test same vs diff patient coverage, P', 'KS_same_vs_diff_pat_coverage', p_ks, '%.1e');

% clamp for log scale
eta0_same  = max(eta0_same,  ETA0_CLAMP);
eta0_other = max(eta0_other, ETA0_CLAMP);
eta0_focal_same = max(eta0_focal_same, ETA0_CLAMP);
eta0_focal_other = max(eta0_focal_other, ETA0_CLAMP);

% focal point positions on clamped CDF (last one after sorting for plot visibility)
ix_s = find(eta0_same  == eta0_focal_same, 1, 'last');
ix_o = find(eta0_other == eta0_focal_other, 1, 'last');

figure(1); clf; hold on;
set(gca, 'xscale', 'log');
% cover the spurious minor tick to the left of ETA0_CLAMP
yl = [0, 1];
ylim(yl);
xlim([ETA0_CLAMP*0.8, 1/0.8]);
plot(ETA0_CLAMP*0.9 + [0 0], yl + [0.002, -0.002], '-w', 'LineWidth', 1);

h_same  = plot(eta0_same,  p_same,  '.-', 'LineWidth', 1, 'Color', csamepat, 'MarkerSize', 9);
h_other = plot(eta0_other, p_other, '.-', 'LineWidth', 1, 'Color', cdiffpat, 'MarkerSize', 9);

plot(eta0_same(ix_s),  p_same(ix_s),  'ok', 'MarkerSize', 4, 'LineWidth', 1);
plot(eta0_other(ix_o), p_other(ix_o), 'ok', 'MarkerSize', 4, 'LineWidth', 1);

% Find x-value of max CDF difference and annotate with significance arrow
all_x = sort([eta0_same; eta0_other]);
cdf_s = mean(eta0_same' <= all_x, 2);
cdf_o = mean(eta0_other' <= all_x, 2);
[~, idx_max] = max(abs(cdf_s - cdf_o));
x_ks = all_x(idx_max);
y_lo = min(cdf_s(idx_max), cdf_o(idx_max));
y_hi = max(cdf_s(idx_max), cdf_o(idx_max));

set(gca, 'XTick', 10.^(-4:0), 'XMinorTick', 'on');

xlabel({'Fraction of urine isolate genome', ...
    'not covered by faecal pathometagenome, {\eta}_0'}, 'FontSize', 7);
ylabel('Cumulative fraction of faecal-urine pairs', 'FontSize', 7);
legend([h_same, h_other], ...
    {'Same-patient', 'Different-patient'}, ...
    'Location', 'nw', 'FontSize', 6, 'Box', 'off');
set(gca, 'FontSize', 6, 'FontName', 'Helvetica', 'Box', 'on');

set(gca, 'Position', [0.2 0.25 0.75 0.7]);
set_axes_physical_size(gca, 2.4, 1.5, 'inch');

h_arrow = annot_arrow(gca, x_ks, y_lo, x_ks, y_hi, 'doublearrow', 'Color', 'k', 'LineWidth', 0.5);
set(h_arrow, 'Head1Length', 4, 'Head1Width', 4, 'Head2Length', 4, 'Head2Width', 4);
text(x_ks, (y_lo + y_hi) / 2, [significance_marker(p_ks) ' '], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 5, 'FontWeight', 'bold');

print_figure(gcf, ['fig_UTI_by_FIT_etha0_' AGGREGATION_LEVEL]);

%% ========================================================================
%  Figure 2: Coverage histogram for focal UTI isolate
%  ========================================================================

BIN_WIDTH = 0.05;
ZERO_VAL = -3;
YMAX = 0.2;

cov_same  = load_and_sum_coverage(runT, output_dir, FOCAL_UTI_ISO, focal_uti.fid);
cov_other = load_and_sum_coverage(runT, output_dir, FOCAL_UTI_ISO, FOCAL_OTHER_FID);
genome_len = numel(cov_same);

cov_same_log  = normalize_and_log(cov_same, ZERO_VAL);
cov_other_log = normalize_and_log(cov_other, ZERO_VAL);

bin_edges = ZERO_VAL : BIN_WIDTH : 1.2;

figure(2); clf; hold on;

h2 = histogram(cov_other_log, 'BinEdges', bin_edges - BIN_WIDTH/2, 'Visible', 'off');
b2 = bar(h2.BinEdges(2:end) - h2.BinWidth*0.5, h2.Values / genome_len, ...
    'FaceColor', cdiffpat, 'EdgeColor', 'none');

h1 = histogram(cov_same_log, 'BinEdges', bin_edges, 'Visible', 'off');
b1 = bar(h1.BinEdges(2:end) - h1.BinWidth*0.5, h1.Values / genome_len, ...
    'FaceColor', csamepat, 'EdgeColor', 'none');

xtick_vals = ZERO_VAL:1;
xtick_labels = arrayfun(@(x) sprintf('10^{%d}', x), xtick_vals, 'UniformOutput', false);
xtick_labels{1} = 'Not covered';

set(gca, ...
    'Box', 'on', ...
    'YLim', [0 YMAX], ...
    'YTick', 0:0.02:YMAX, ...
    'XLim', [ZERO_VAL - 0.1, 0.5], ...
    'XTick', xtick_vals, ...
    'XTickLabel', xtick_labels, ...
    'XTickLabelRotation', 45, ...
    'FontName', 'Helvetica', ...
    'FontSize', 6);

uti_isolate = char(focal_uti.sample_name);
ylabel({'Fraction of urine isolate genome', 'covered by faecal pathometagenome'});
xlabel({['Coverage of urine isolate ', uti_isolate], '(normalized to median)'});

same_label  = sprintf('Same patient (P%d)',  focal_uti.fid);
other_label = sprintf('Different patient (P%d)', FOCAL_OTHER_FID);
legend([b1, b2], {same_label, other_label}, ...
    'Box', 'off', 'Location', 'nw', 'FontSize', 6);

set(gca, 'Position', [0.15 0.3 0.8 0.65]);
set_axes_physical_size(gca, 2.4, 2.0, 'inch');
print_figure(gcf, 'fig_UTI_by_FIT_coverage_histogram');

end


%% ---- local helpers -----------------------------------------------------

function cov_sum = load_and_sum_coverage(runT, output_dir, uti_iso, fid)
    mask = runT.uti_run_name == uti_iso & ...
           runT.faecal_fid == fid & ...
           runT.status == "ok";
    rows = runT(mask, :);
    assert(height(rows) >= 1, 'No coverage runs found for %s / F%d', uti_iso, fid);

    cov_sum = [];
    for i = 1:height(rows)
        npz_path = fullfile(output_dir, 'alignment', rows.run_name(i), 'coverage.npz');
        cov = load_npz_coverage(char(npz_path));
        if isempty(cov_sum)
            cov_sum = cov;
        else
            cov_sum = cov_sum + cov;
        end
    end
end


function cov = load_npz_coverage(npz_path)
    tmp = [tempname '.bin'];
    cmd = sprintf(['python3 -c "import numpy as np; ' ...
        'c=np.load(''%s'')[''coverage'']; ' ...
        'c.astype(''<i4'').tofile(''%s'')"'], npz_path, tmp);
    [status, msg] = system(cmd);
    assert(status == 0, 'Failed to load %s: %s', npz_path, msg);
    fid = fopen(tmp, 'r');
    cov = fread(fid, Inf, 'int32');
    cov = double(cov);
    assert(all(cov >= 0));
    fclose(fid);
    delete(tmp);
end


function norm_cov = normalize_and_log(cov, zero_val)
    med = median(cov);
    norm_cov = log10(cov / med);
    norm_cov = max(norm_cov, zero_val);
end


function T = read_stats_table(output_dir, filename)
    T = readtable(fullfile(output_dir, 'summary', filename), 'TextType', 'string');
end
