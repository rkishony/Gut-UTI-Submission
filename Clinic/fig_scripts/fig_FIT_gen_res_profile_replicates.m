function [cf, h] = fig_FIT_gen_res_profile_replicates(fig, datasets, p_vals)

cf = [];
h = [];

global CONFIG  %#ok<GVMIS>
if ~CONFIG.create_figures
    return
end

rng(43)

clr = [0.5 0.5 0.5];
jitter_x = 0.05;
jitter_y = 0.02;

grp = repelem(1:3, cellfun(@numel, datasets));

cf = figure(fig);
clf
hold on; box on;

h.bosplot = boxplot([datasets{:}], grp, 'colors',clr, 'symbol','');
axis([0.5 3.5 -0.1 1.1])
for idata = 1:3
    dataset = datasets{idata};
    
    if idata~=3
        jit_x = jitter_x * randn(numel(dataset),1);
        jit_y = jitter_y * rand(numel(dataset),1);
        h.scatter = scatter(idata + jit_x, dataset' + jit_y, 14, ...
            'marker','o','markerfacecolor','k','markeredgecolor','none',...
            'markerfacealpha',0.5);
    end

    % h.line = line(idata + jitter_x*[-1 1], mean(dataset) + [0 0], 'color','k','linewidth',2);

    h.n_eq_text = text(idata, 1.1, sprintf('n=%d',length(datasets{idata})), ...
        'fontsize',6,'HorizontalAlignment','center','VerticalAlignment','bottom');
end
set(gca,'fontsize',7,'position',[0.15 0.2 0.8 0.75],'xtick',[1 2 3],...
    'TickLabelInterpreter','tex',...
    'XTickLabel',{'Same-patient,\newlineReplicate cohort','Cross-patient,\newlineReplicate cohort','Cross-patient,\newlineAll'})
ylabel('Distance between genotype-inferred resistance profiles (Hamming)');

sigstar({[1 2], [1 3], [2 3]}, [p_vals.p12, p_vals.p13, p_vals.p23])

set_axes_physical_size(gca,4.5,3.5,'inch',1);
end
