function ld  = fig_legend(fig, varargin)
% like legend, but can plot in a new figure if fig is provided

ld = legend(varargin{:});
if isempty(fig)
    return
end

set(ld, 'units', 'norm')
ori_ld_norm_pos = get(ld,'Position');

ori_fig = ancestor(ld, 'figure');
flds = {'PaperUnits', 'PaperPosition', 'Units', 'Position'};
vals = get(ori_fig, flds);
fig_props = cell2struct(vals, flds, 2);


fig = figure(fig);
set(fig,'ToolBar','none', 'MenuBar','none')
clf

ax = gca;
set(ax,'visible','off')
ld = legend(ax, varargin{:});
fig_props.PaperPosition(3:4) = fig_props.PaperPosition(3:4) .* ori_ld_norm_pos(3:4);
fig_props.Position(3:4) = fig_props.Position(3:4) .* ori_ld_norm_pos(3:4);
set(fig, fig_props);
set(ld, 'units', fig_props.PaperUnits)
set(ld, 'position', [0 0 fig_props.PaperPosition(3:4)])

end