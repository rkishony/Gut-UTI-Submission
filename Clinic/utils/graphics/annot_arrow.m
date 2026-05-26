function h = annot_arrow(ax, x1, y1, x2, y2, type, varargin)
% Adds an annotation arrow from (x1,y1) to (x2,y2) in axis data units,
% supports linear and log scales.

if nargin < 6
    type = 'arrow';
end

fig = ancestor(ax, 'figure');
ax_pos = get(ax, 'Position');  % [left bottom width height] in normalized units

% Get axis limits
ax_xlim = xlim(ax);
ax_ylim = ylim(ax);

% Handle log scales
if strcmp(get(ax, 'XScale'), 'log')
    x1 = log10(x1);
    x2 = log10(x2);
    ax_xlim = log10(ax_xlim);
end
if strcmp(get(ax, 'YScale'), 'log')
    y1 = log10(y1);
    y2 = log10(y2);
    ax_ylim = log10(ax_ylim);
end

% Normalize
x_norm = ([x1, x2] - ax_xlim(1)) / diff(ax_xlim);
y_norm = ([y1, y2] - ax_ylim(1)) / diff(ax_ylim);

% Convert to figure units
x_fig = ax_pos(1) + ax_pos(3) * x_norm;
y_fig = ax_pos(2) + ax_pos(4) * y_norm;

% Create annotation
h = annotation(fig, type, x_fig, y_fig, varargin{:});
end
