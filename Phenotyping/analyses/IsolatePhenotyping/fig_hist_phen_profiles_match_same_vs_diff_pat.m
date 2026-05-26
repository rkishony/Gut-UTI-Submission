function h = fig_hist_phen_profiles_match_same_vs_diff_pat(fignum, b_same, b_bs, nbins)

if nargin<4
    nbins = 36;
end
figure(fignum); clf; hold on
set(gcf,'PaperPosition',[0 0 3 3], 'Units','inches')
set(gca,'fontsize',5,'FontName','Helvetica')

permutedclr = [0.65 0.65 0.65];

h.hist = histogram(b_bs(1,:), linspace(0,0.6,nbins), 'FaceColor',permutedclr , 'FaceAlpha',1, 'EdgeColor',[1 1 1], 'EdgeAlpha',0.5);%, 'Orientation', 'horizontal');
arrowtop = max(h.hist.BinCounts)/5;
h.yl = plot(b_same(1)+[0 0], [arrowtop 48],'linewidth',4,'Color','#e2813e');
h.ar = plot(b_same(1), 48, 'v','MarkerFaceColor','#e2813e', 'MarkerEdgeColor', '#e2813e', 'MarkerSize',8);
xlabel('Fraction of pathobiome-urine sample pairs with a perfect match');
ylabel('Number of permutations');
xlim([0, 0.55])
box on

h.text(1) = text(b_same(1), arrowtop*1.05, 'Same patient','HorizontalAlignment','right','VerticalAlignment','baseline','FontSize',5);
% h.text(1) = text(0.6*0.999, arrowtop*1.1, 'Same patient','HorizontalAlignment','right','VerticalAlignment','baseline','FontSize',5);
h.text(2) = text(mean(b_bs(1,:)), max(h.hist.BinCounts)+20, 'Permuted', 'HorizontalAlignment','center','VerticalAlignment','baseline','FontSize',5);

% set(gca, 'YLim', [0 2100])
set(gca, 'YLim', [0 round(max(h.hist.Values)*1.1,-1)])

n_rand = length(b_bs(1,:));
pval = length(find(b_bs(1,:)>=b_same(1)))/n_rand;
if pval==0
    txt_pval = sprintf('p<%.4f', 1/n_rand);
else
    txt_pval = sprintf('p=%.4f', pval);
end
h.text(3) = text(0.55*0.995, max(get(gca,'YLim'))*0.995, txt_pval, 'fontsize', 5, 'HorizontalAlignment','right','VerticalAlignment','top');

end
