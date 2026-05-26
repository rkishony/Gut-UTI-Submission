function sz = get_physical_size(ax)
fig=ax.Parent;
axpos=get(ax,'Position');
figpos=get(fig,'PaperPosition');
sz=[axpos(3:4).*figpos(3:4)];
end
