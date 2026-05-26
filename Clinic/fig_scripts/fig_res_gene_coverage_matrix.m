function [cf, num_samples_not_plotted] = fig_res_gene_coverage_matrix(fig, RF, RF_Entities, ent_name, grp_name, thresh)

cf = [];

% Coverage 
cov = RF.(['cov_' ent_name]);

% Choose samples (remove samples with no coverage on any genes)
jSamples = find(sum(cov>thresh,2) > 0);
nSamples = numel(jSamples);
toplot = cov(jSamples,:);
num_samples_not_plotted = size(cov,1) - size(toplot,1);

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end

% Remove genes with no coverage in any samples
jEntity = find(sum(toplot>thresh,1)>0);
nEntity = numel(jEntity);
toplot = toplot(:,jEntity);

% Sort the samples by clustering similarity:
im = sort_rows_by_cluster(double(toplot>thresh));
toplot = toplot(im,:);

% Anything smaller than threshold should appear as "0" (white):
% toplot(toplot<thresh) = 0;
toplot = log10(toplot);

% Figure
cf = figure(fig); clf
ylabel('Resistance genes','FontSize',7);
xlabel('Faecal pathometagenomes','FontSize',7);

% Colormap
c = linspace(0.95, 0, 40);
cmap = c' + [0 0 0];

% Plot the image
h = imagesc(toplot');
hold on
h.AlphaData = toplot' > log10(thresh);
clim([log10(thresh), log10(8)])
set(gca,'xTick',1:nSamples,'ticklength',[0.002 0],'xticklabelrotation',90,...
    'xticklabel',RF.Label(jSamples(im)),...
    'YTick',1:nEntity,...
    'YTickLabel',RF_Entities.(['cov_' ent_name]){jEntity, ent_name},...
    'fontsize',3.5, 'colormap',cmap, ...
    'Position',[0.05 0.33 0.85 0.65]);
xlabel('Patient faecal pathobiome', 'fontsize',7)
ylabel('Resistance genes', 'fontsize',7)

% Add seperator lines and group names:
if ~isempty(grp_name)
    grp_ent = RF_Entities.(['cov_' ent_name]){jEntity, grp_name};
    grpend = find(~strcmp(grp_ent(2:end), grp_ent(1:end-1)));
    yline(grpend+0.5,'k-')
    grpend = [0; grpend; numel(jEntity)];
    for i = 1:numel(grpend)-1
        nm = grp_ent(grpend(i+1));
        if strcmp(nm, 'beta_lactam')
            nm = '\beta-lactam';
        end
        x0 = nSamples+1;
        y0 = (grpend(i) + grpend(i+1))/2+0.5;
        x1 = x0;
        y1 = y0;
        if strcmp(nm, 'sulphonamide')
            x1 = x1 + 1.3;
            y1 = y1 - 0;
        end
        if strcmp(nm, 'tetracycline')
            x1 = x1 + 1.3;
            y1 = y1 + 0.7;
        end
        if strcmp(nm, 'trimethoprim')
            x1 = x1 + 1.3;
            y1 = y1 + 2 * 0.7;
        end

        text(x1, y1, nm, 'horizontalalignment','left','verticalalignment','middle', 'fontsize',6)
        if x1 > x0
            plot([x0-0.02, x1-0.2], [y0 y1], 'k-', 'linewidth', 0.5, 'clipping','off')
        end
    end
end

% Colorbar
cb=colorbar;
ticks = [0.05 0.1 0.2 0.4 0.8 1 2 4 8];
sticks = arrayfun(@num2str, ticks, 'uni', 0);
set(cb,'xtick',log10(ticks), 'xticklabel',sticks, ...
    'Position',[0.6 0.15 0.3 0.02], 'FontSize',5, 'Orientation','horizontal')
ylabel(cb, 'Resistance gene coverage, normalized','fontsize',7);

set_axes_physical_size(gca,6,6/diff(xlim)*diff(ylim)*1.25,'inch');

% "<" white box
ax = axes('Position',[0.57, 0.15, 0.008, 0.02], 'FontSize',5);
axis([0 1 0 1])
set(gca,'xtick',[],'ytick',[])
set(ax,'box','on')
text(0.8, -0.5, sprintf('<%g',thresh), 'fontsize', 5, 'clipping', 'off', ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle', 'Rotation',60)


end