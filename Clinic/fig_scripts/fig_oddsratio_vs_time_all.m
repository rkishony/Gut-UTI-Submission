function cf = fig_oddsratio_vs_time_all(fig, OR, OR_perm_P_val, OR_perm_prctile, OR_bs_prctile, ...
    Days, legend_labels)

cf = [];

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end

if nargin<8
    plot_perm_null = false;
end

inf_val = log2(10);

% Convert to log:
OR = bounded_log2(OR, 2^inf_val);
OR_perm_prctile = bounded_log2(OR_perm_prctile, 2^inf_val);
OR_bs_prctile = bounded_log2(OR_bs_prctile, 2^inf_val);

order = [1 3 2];
OR = permute(OR, order);
OR_perm_prctile = permute(OR_perm_prctile, order);
OR_bs_prctile = permute(OR_bs_prctile, order);
OR_perm_P_val = permute(OR_perm_P_val, order);

% Define Confidence Interval to plot:
OR_bs_L = OR_bs_prctile(:,:,1);
OR_bs_H = OR_bs_prctile(:,:,3);

% Time axis
Days = Days';
Weeks = Days / 7;

% Figure
cf = figure(fig);clf
hold on; box on
set(gca, 'Position', [0.1 0.08 0.88 0.7])
xlabel('Time post faecal sample, weeks', 'FontSize',8)
ylabel({'Odds ratio for urine culture resistance', ...
    'given faecal resistance'}, 'FontSize',8)
xlim([-4 44])
ylim([0 log2(11)])

set(gca,'XTick',Weeks)
ytck = 1:10;
ytl = arrayfun(@num2str, ytck, 'uni', 0);
ytl{end} = ['>' ytl{end}];
set(gca,'Ytick',log2(ytck),'YtickLabel',ytl)
set(gca,'FontSize',7)

% Nominal values:
hbar = bar(Weeks, OR, 'BarWidth',0.8);
x = arrayfun(@(h) h.XEndPoints, hbar, 'UniformOutput',false);
x = cat(1, x{:})';
for h = hbar
    h.FaceColor = 1 - (1 - h.FaceColor)*0.3;
end

% Error bars:
OR_bs_H(OR_bs_L<0) = OR(OR_bs_L<0);
OR_bs_L(OR_bs_L<0) = OR(OR_bs_L<0);
errorbar(x, OR, OR-OR_bs_L, -OR+OR_bs_H, 'Color','k', ...
    'LineStyle','none')

% Significance:
for jt = 1:numel(OR)
    y = OR_bs_H(jt) + 0.02;
    text(x(jt), y, significance_marker(OR_perm_P_val(jt)), 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',5, 'Clipping','on')
end

set_axes_physical_size(gca, 5.5, 3.5, 'inches');

% legend
normal_axes_legend(cf, [0.1 0.8 0.3 0.15], hbar, legend_labels, 'fontsize',6);


end


function log2x = bounded_log2(x, maxval)
if nargin<2
    maxval = 1e8;
end
x(x>maxval) = maxval;
x(x<1/maxval) = 1/maxval;
log2x = log2(x);
end