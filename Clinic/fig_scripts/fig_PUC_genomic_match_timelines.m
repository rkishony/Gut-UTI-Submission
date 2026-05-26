function cf = fig_PUC_genomic_match_timelines(fig, UTIpop, SMP, samepatcovbreadth, etha0_th)

cf = [];

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end

samepatcovbreadth = addMatch(samepatcovbreadth, UTIpop, 'UTI_ID', 'UTI_Drisha');
samepatcovbreadth = addMatch(samepatcovbreadth, SMP, 'UTI_Drisha', {'FOBT_ID', 'DateDiff'});
coveredutis = samepatcovbreadth.etha0 < etha0_th;

[u_ids,~,jid] = unique(samepatcovbreadth.FOBT_ID);
max_t = max(samepatcovbreadth.DateDiff) + 20;

clrs = [
    0.5 0.5 0.5
    0.8 0.6 0
    ];


cf = figure(fig);
clf
hold on; box on;
set(gca,'fontsize',6)
axis([-9 max_t+10 0 length(u_ids)+1]);

for i = 1:length(u_ids)
    plot([-7 max_t],[i i],'k-','LineWidth',0.5);
    plot(max_t, i, 'k>', markersize=6, MarkerFaceColor='k');
end
set(gca,'ytick', 1:length(u_ids),'yticklabel',arrayfun(@num2str, u_ids, 'UniformOutput', false));
set(gca,'xtick', 0:14:200, 'xticklabels', arrayfun(@num2str, (0:14:200)/7, 'UniformOutput', false))
ylabel('Patient ID', 'FontSize',7)
xlabel('Time since faecal FIT sample, weeks', 'fontsize',7)


for i = 1:height(samepatcovbreadth)
    plot(samepatcovbreadth.DateDiff(i), jid(i)+0.5, 'v', 'MarkerSize',7, 'markeredgecolor','none',...
        'MarkerFaceColor', clrs(coveredutis(i)+1,:))
end

set_axes_physical_size(gca,4.5,5,'inch',1);

end