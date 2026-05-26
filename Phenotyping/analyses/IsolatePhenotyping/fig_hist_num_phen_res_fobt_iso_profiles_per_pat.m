function h = fig_hist_num_phen_res_fobt_iso_profiles_per_pat(fignum, nprofs)
figure(fignum); clf; 
hold on;

edges = (min(nprofs)-0.5):1:(max(nprofs)+0.5);  % bin edges for integer bins
counts = histcounts(nprofs, edges);

binCenters = edges(1:end-1) + diff(edges)/2;
blockGap = 0.03;        % vertical spacing between units
blockHeight = 1-blockGap;
barWidth = 0.8;         % width of each block
yShift = 0;             % optional y offset

for i = 1:length(counts)
    x = binCenters(i);
    for j = 1:counts(i)
        yBottom = (j-1)*(blockHeight + blockGap) + yShift;
        h.rects(i,j) = rectangle('Position', [x - barWidth/2, yBottom, barWidth, blockHeight], ...
                  'FaceColor', [0.5 0.5 0.5], ...
                  'EdgeColor', [0.75 0.75 0.75]);
    end
end

xlim([edges(1)-0.25, edges(end)+0.25]);
ylim([0 round(max(counts)*1.05)]);
xlabel('Number of unique profiles');
ylabel('Number of patients');
box on
set(gcf,'PaperPosition',[0 0 4.5 2.5],'Units','inches')
set(gca,'fontsize',5,'FontName','Helvetica')

% make highest rectangle on fourth column red
highestfourth = find(isgraphics(h.rects(4,:)),1,'last');
set(h.rects(4,highestfourth), 'FaceColor', '#91171f')

end