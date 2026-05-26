function [p_values, X, Y, h_lines, h_arrows, h_texts] = plot_cum_hists(ax, chosen, xs, vertical, x_c, varargin)

if nargin<4
    vertical = 'doublearrow';
end
if nargin<5
    x_c = [];
end

is_log = strcmpi(get(gca,'xscale'),'log');
arrow_head = 4;

xs = xs(:);
hold(ax, 'on');
n_groups = size(chosen, 2);
assert(size(chosen,1)==size(xs,1));
assert(iscolumn(xs))
h_lines = gobjects(n_groups,1);
h_arrows = gobjects(n_groups-1,1);
h_texts = gobjects(n_groups-1,1);
p_values = nan(n_groups-1,1);

% Expand style args into per-group args
per_group_args = expand_style_args(varargin, n_groups);
is_log = strcmp(get(gca,'xscale'),'log');
x_lim = xlim;

X = nan(size(chosen));
Y = nan(size(chosen));

for i = 1:n_groups
    csn = chosen(:,i);
    x_i = xs(csn);
    [P_i, x_i, k] = getP(x_i);
    X(csn, i) = x_i(k);
    Y(csn, i) = P_i(k);
    h_lines(i) = plot(ax, x_i, P_i, per_group_args{i}{:});
    if i==1
        x_r = x_i; 
        P_r = P_i;
    else
        [~, p_val] = kstest2(x_r, x_i);
        p_values(i-1) = p_val;

        % max dist location
        x_all = sort([x_r; x_i]);
        P_all_ref = mean(x_r' <= x_all, 2);
        P_all_i = mean(x_i' <= x_all, 2);
        [~, idx] = max(abs(P_all_ref - P_all_i));
        x_at_max = x_all(idx);
        if isempty(x_c)
            x_arrow = x_at_max;
        else
            shift = - (n_groups - i + 1) * 0.008;
            if is_log
                x_arrow = exp(log(x_c) + diff(log(x_lim)) * shift);
            else
                x_arrow = x_c + diff(x_lim) * shift;
            end
        end
        y1 = cdf_interp(x_r, P_r, x_arrow, is_log);
        y2 = cdf_interp(x_i, P_i, x_arrow, is_log);
        if ~strcmp(vertical,'none')
            switch vertical
                case 'line'
                    h_arrows(i-1) = plot(ax, [x_arrow x_arrow], [y1 y2], 'k--x');
                case {'arrow', 'doublearrow'}
                    h_arrow = annot_arrow(ax, x_arrow, y1, x_arrow, y2, vertical, 'color', 'k', 'LineWidth', 1);
                    if strcmp(vertical,'arrow')
                        set(h_arrow, 'HeadLength', arrow_head, 'HeadWidth', arrow_head)
                    else
                        set(h_arrow, 'Head1Length', arrow_head, 'Head1Width', arrow_head)
                        set(h_arrow, 'Head2Length', arrow_head, 'Head2Width', arrow_head)
                    end
                    set(h_arrow, 'linewidth', 0.5)
                    h_arrows(i-1) = h_arrow;
            end
            h_texts(i-1) = text(ax, x_arrow, y2+0.02, [significance_marker(p_val) ' '], ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 5, 'FontWeight', 'bold');
        end
    end
end

end

function argsets = expand_style_args(args, n)
argsets = cell(n,1);
for i = 1:n
    arg_i = {};
    for j = 1:2:length(args)
        key = args{j};
        val = args{j+1};
        if iscell(val)
            val_i = val{i};
        else
            val_i = val;
        end
        arg_i = [arg_i, {key, val_i}];
    end
    argsets{i} = arg_i;
end
end


function [P, x, ind] = getP(x)
n = numel(x);
[x, ind] = sort(x); 
[~,ind ] = sort(ind);
% P = (1:n)/n;
P = linspace(0, 1, n);
end

function y0 = cdf_interp(x,y,x0,is_log)
dx = diff(x);
dx(dx==0) = 1e-8;
x = x(1)+[0; cumsum(dx)];
%[x, i] = unique(x);
%y = y(i);
if is_log
    y0 = interp1(log(x),y,log(x0),'linear',0);
else
    y0 = interp1(x,y,x0,'linear',0);
end

end
