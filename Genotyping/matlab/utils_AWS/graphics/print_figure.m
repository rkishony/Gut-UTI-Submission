function print_figure(fig, name, fmts, correct_svg_size, enforce_font)

global CONFIG %#ok<GVMIS>

if nargin<5
    enforce_font = [];
end
drawnow
if ~isempty(fig)
    set(fig, 'Name', name)
    if ~isempty(enforce_font)
        set(findall(gcf, '-property', 'FontName'), 'FontName', enforce_font);
    end
else
    return
end

should_print = true;
if ~isempty(CONFIG) && isstruct(CONFIG) && isfield(CONFIG, 'print_figures')
    should_print = logical(CONFIG.print_figures);
end
if ~should_print
    return
end

if nargin<3 || isempty(fmts)
    fmts = 'svg';
end
if nargin<4 || isempty(correct_svg_size)
    correct_svg_size = true;
end

if ~iscell(fmts)
    fmts = {fmts};
end

folders = get_folders();
if ~exist(folders.figures, 'dir')
    mkdir(folders.figures);
end
filepath = fullfile(folders.figures, name);

clr = get(fig,'Color');
set(fig,'Color','w')
set(fig, 'InvertHardcopy', 'off')
set(findall(gcf,'-property','FontName'),'FontName','Arial');

for i = 1:numel(fmts)
    fmt = fmts{i};
    print(filepath, '-vector', ['-d' fmt])
    if strcmp(fmt, 'svg') && correct_svg_size
        old_units = fig.PaperUnits;
        fig.PaperUnits = 'inches';
        pos = fig.PaperPosition;
        fig.PaperUnits = old_units;
        fix_svg_size([filepath '.svg'], pos(3), pos(4))
    end
end
set(fig,'Color',clr)

end

