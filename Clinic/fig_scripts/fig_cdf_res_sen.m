function [cf, h, X, Y] = fig_cdf_res_sen(...
    fig, rf_or_phen, fobt_val, chosen, pretty, fobt_val_thr, is_log, x_lim, clr_sen, clr_res)

cf = [];
h = [];
X = [];
Y = [];

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end


assert(iscolumn(fobt_val))

% Set labels
legends = {
    'Only sensitive urine cultures'
    '>= 1 resistant urine cultures'
    '>= 2 resistant urine cultures'
    '>= 3 resistant urine cultures'
    };

if strcmp(rf_or_phen,'RF')
    if pretty == "CPR (mut)"
        xlbl = {sprintf('Number of CPR resistant mutations in', pretty), 'patient''s pathometagenome, \itf_{genotype}'};
    else
        xlbl = {sprintf('Normalized coverage of %s genes in', pretty), 'patient''s pathometagenome, \itf_{genotype}'};
    end
    fc_txt = ' \itf_{genotype}^{ C}';
else
    xlbl = {sprintf('Frequency of %s resistant colonies in', pretty), 'patient''s pathobiome, \itf_{phenotype}'};
    fc_txt = ' \itf_{phenotype}^{ C}';
end

% Set line styles
shadeconstant = 0.7; 
colors = num2cell([clr_sen; clr_res; clr_res; clr_res]*shadeconstant,2);
markers = num2cell('..*+');
markersizes = {8, 8, 4, 4};
n = size(chosen,2);
colors = colors(1:n);
markers = markers(1:n);
markersizes = markersizes(1:n);
legends = legends(1:n);

% Set x-lim
if isnan(x_lim(1))
    x_lim(1) = min(fobt_val);
end
if isnan(x_lim(2))
    x_lim(2) = max(fobt_val);
end

fobt_val(fobt_val<x_lim(1)) = x_lim(1);
fobt_val(fobt_val>x_lim(2)) = x_lim(2);

% extend x_lim a bit:
if is_log
    dx = diff(log(x_lim));
    x_lim_actual = exp(log(x_lim) + dx*0.01*[-1 1]);
else
    dx = diff(x_lim);
    x_lim_actual = x_lim + dx*0.01*[-1 1];
end

% Create figure
cf = figure(fig);
clf
hold on
ylim([0 1])
xlim(x_lim_actual)
if is_log
    set(gca, 'xscale','log')
    xticks = floor(log10(x_lim(1))):ceil(log10(x_lim(2)));
    set(gca, 'xtick',10.^xticks)
end
set(gca, 'box','on', 'Position',[0.15, 0.25, 0.8, 0.67], 'fontsize',6)

[KS_p_values, X, Y, h.lines, h.arrows, h.texts] = plot_cum_hists(gca, chosen, fobt_val, 'arrow', fobt_val_thr, ...
    'color',colors, 'marker',markers, 'markersize',markersizes);
xlabel(xlbl,'interpreter','tex', 'fontsize',7)
ylabel('Cumulative fraction of patients', 'fontsize',7);

% Threshold line
h.val_th_line = xline(fobt_val_thr,'color','black','linewidth',0.5,'linestyle','--');
h.val_th_text = text(fobt_val_thr, 0.06, fc_txt, 'HorizontalAlignment','left', 'VerticalAlignment','middle','FontSize',6);

h.rectangle(1) = rectangle('Position',[x_lim(1) 1.02 fobt_val_thr-x_lim(1) 0.05],'clipping','off','edgecolor','none','facecolor',clr_sen);
h.rectangle(2) = rectangle('Position',[fobt_val_thr 1.02 x_lim(2)-fobt_val_thr 0.05],'clipping','off','edgecolor','none','facecolor',clr_res);

h.legend = legend([plot(nan,nan,'-','color','none'); h.lines], ...
    [{'Patients with:'}; legends], 'fontsize', 6, ...
    'Position',[0.26 0.75 0.3 0.15], 'box','off');
h.legend.AutoUpdate = "off";

if rf_or_phen == "PHEN"
    set(legend, 'Location', 'southeast')
end

set_axes_physical_size(gca, 3.2, 2.5, 'inch', 2);

end
