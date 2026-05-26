function [figW, figH] = set_axes_physical_size(ax, target_width, target_height, units, screen_factor)
if nargin<4, units = 'centimeters'; end
if nargin<5, screen_factor = 2; end

fig = ancestor(ax, 'figure');

% assert(strcmp(get(gca,'Units'), 'norm')
set(fig, 'PaperUnits', units);
set(fig, 'Units', units)

pos = get(ax, 'Position');

figW = target_width / pos(3);
figH = target_height / pos(4);

set(fig, 'PaperPosition', [0 0 figW, figH])
set(fig, 'Position', [1, 1, figW, figH] * screen_factor)
warning off
set(ax, 'Position', pos)
warning on
end
