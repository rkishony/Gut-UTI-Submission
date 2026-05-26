function [X, h_bar, h_text] = multidimbar(Y, ws, labels, lbl_offset)
if nargin<3
    labels = [];
end
if nargin<4
    lbl_offset = 0.04;
end

szs = size(Y);
ndim = numel(szs);
w1 = ws(1);

ctrs = cell(size(szs));
pw = 1;
for d = ndim:-1:1
    sz = szs(d); w = ws(d);
    gap = (pw - sz*w) / sz;
    ctrs{d} = (w+gap) * ((1:sz) - sz/2 - 0.5);
    pw = w;
end

X = zeros(szs);
h_bar = gobjects(szs);

subsMat = zeros(numel(Y), ndim);
for i = 1:numel(Y)
    subs = ind2subvec(szs,i);
    subsMat(i,:) = subs;

    ctr = 0.5;
    for d = 1:ndim
        ctr = ctr + ctrs{d}(subs(d));
    end
    X(i) = ctr;
    h_bar(i) = rectangle( ...
        'Position',[ctr-w1/2, 0, w1, Y(i)], ...
        'FaceColor','b', 'EdgeColor','none');

end
xlim([0 1])

% Add x-labels
yl = ylim();
labelOffset = lbl_offset * diff(yl);
set(gca,'xtick',sort(X(:)),'XTickLabel',[])

level = 0;
h_text = cell(1,ndim);
for d = 1:ndim
    if numel(labels)<d || isempty(labels{d})
        h_text{d} = gobjects([1 0]);
        continue
    end

    dimsIdx = d:ndim;
    nC      = numel(dimsIdx);

    gridOut = cell(1,nC);
    C = arrayfun(@(k) 1:szs(dimsIdx(k)), 1:nC, 'uni',false);
    [gridOut{:}] = ndgrid(C{:});

    tmp   = cellfun(@(g) g(:), gridOut, 'uni',false);
    combs = cat(2, tmp{:});
    hs = gobjects(1, size(combs,1));
    for r = 1:size(combs,1)
        sel  = all(bsxfun(@eq, subsMat(:,dimsIdx), combs(r,:)),2);
        xMed = mean(X(sel));
        yPos = -(level+0.1) * labelOffset;
        txt  = labels{d}{combs(r,1)};
        hs(r) = text(xMed, yPos, txt, 'HorizontalAlignment','center','VerticalAlignment','top');
    end
    h_text{d} = hs;
    level = level + 1;
end

end


function subs = ind2subvec(sz, idx)
n = numel(sz);
[subs_cell{1:n}] = ind2sub(sz, idx);
subs = [subs_cell{:}];
end