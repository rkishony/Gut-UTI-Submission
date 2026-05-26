%% plot IsolatePhenMatch vs CovBreadth
function plot_IsolatePhenMatch_vs_CovBreadth(hbin, binedges, plotname)
lg = {'0 (perfect match)','1','>=2'};
fN = figure('Name','Histogram of Analyzed Breadth to Phenotype'); hold on
xdata = repmat((binedges(1:end-1)+binedges(2:end))/2,3,1)';
ydata = hbin;
b = bar(xdata, ydata, 'stacked', 'BarWidth', 1);

clr_identical = '#e48336';
clr_over1 = '#81532c';
clr_diff1 = '#a56026';

b(1).FaceColor = clr_identical; b(1).EdgeColor = clr_identical;
b(2).FaceColor = clr_diff1; b(2).EdgeColor = clr_diff1;
b(3).FaceColor = clr_over1; b(3).EdgeColor = clr_over1;

ylim([0 max(max(b(3).YEndPoints),max(b(3).YEndPoints))+2])
fN.Position = [500   400   850   630];
leg = legend(lg,'location','northwest','Box','off');

leg.Title.String = 'Distance between resistance profiles';
xt = [-4:0];
xticks(xt);
xlim([-4.5 0.5]);
xtl = arrayfun(@(t) sprintf('10^{%d}',t), xt, 'UniformOutput', false);
xtl{1} = ['<=' xtl{1} ' '];
xticklabels(xtl)
% xtickangle(45)
set(gcf, 'color','white');
xl = xlabel({'Fraction of urine isolate genome not covered by','pathometagenome, \eta_0'});

ylabel([{'Number of same-patient', 'pathobiome-urine sample pairs'}])
ax = gca;
ax.FontSize = 6;
leg.Position = [0.23   0.7442    0.3671    0.1618];
box on
set_axes_physical_size(gca, 3, 2.2, 'inch',2);
print_figure(gcf, plotname)
