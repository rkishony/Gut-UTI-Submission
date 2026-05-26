function cf = fig_pop_stat(fig, POP, PUC, isCoinciding, isNonCoinciding)

cf = [];

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end

cf = figure(fig); clf
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact')

% Age
edges = [0 30:10:90 inf];
lbls = arrayfun(@(a,b) sprintf('%d-%d',a,b),edges(1:end-1),edges(2:end),'UniformOutput',false);
lbls{1} = ['<' num2str(edges(2))];
lbls{end} = ['\geq' num2str(edges(end-1))];
plot_two_hists(POP.age, isCoinciding, isNonCoinciding, edges);
xlabel('Age')
ylabel('Fraction of patients')
set(gca,'xtick',1:length(edges)-1,'XTickLabel',lbls)

% Sex
edges = 0:2;
plot_two_hists(POP.Gender=='F', isCoinciding, isNonCoinciding, edges);
xlabel('Sex')
ylabel('Fraction of patients')
set(gca,'xtick',1:2, 'XTickLabel',{'Male', 'Female'})
legend({'Coinciding','Not coinciding'},'Location','nw','box','off')

% Num of PUC
edges = [0:7, inf];
lbls = arrayfun(@num2str,edges(1:end-1),'UniformOutput',false);
lbls{end} = ['\geq' lbls{end}];
num_puc_per_patient = countMatches(POP, PUC, 'FOBT_ID');
plot_two_hists(num_puc_per_patient, isCoinciding, isNonCoinciding, edges);
set(gca,'xtick',1:numel(edges)-1,'XTickLabel',lbls)
xlabel('Number of PUCs per patient')
ylabel('Fraction of patients')

% Date diff PUC - FOBT
edges = [-1 7 28:28:336 366];
[~,jPOP] = ismember(PUC.FOBT_ID, POP.FOBT_ID);
plot_two_hists(PUC.DateDiff, isCoinciding(jPOP), isNonCoinciding(jPOP), edges);
xlim([0.5 length(edges)-0.5])
set(gca,'xtick',0.5:length(edges)-0.5,'XTickLabel',[0, num2cell(edges(2:end-1)/7), '1 year'])

xlabel('Date difference between FIT and PUC (weeks)')
ylabel('Proportion of PUCs')

set_axes_physical_size(gca, 2.8, 1.8, 'inch', 1);

end

function [ax, hbar] = plot_two_hists(val, ind1, ind2, edges)
ax = nexttile;
hold on; box on;
set(gca,'fontsize',6);
c1 = histcounts(val(ind1), edges);
c2 = histcounts(val(ind2), edges);
hbar = bar([c1'/sum(c1), c2'/sum(c2)]);
hbar(1).FaceColor = [0.5 0.5 0.5];
hbar(2).FaceColor = [0.8 0.1 0.1];
end
