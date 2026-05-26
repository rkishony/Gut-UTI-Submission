function rgb = get_colors(names, default_rgb)
%GET_COLORS Return RGB colors (0-1) from project colors.txt by name.
%   rgb = GET_COLORS(names) returns an N-by-3 matrix for the provided
%   names (string/cellstr/char). Missing names are assigned default gray.
%
%   rgb = GET_COLORS(names, default_rgb) sets the fallback color.

if nargin < 2 || isempty(default_rgb)
    default_rgb = [200, 200, 200] ./ 255;
end

if ischar(names)
    names = {names};
elseif isstring(names)
    names = cellstr(names);
end

this_file = mfilename('fullpath');
project_root = fullfile(fileparts(fileparts(this_file)), '..');
colors_file = fullfile(project_root, 'colors.txt');

T = readtable(colors_file, 'Delimiter', ',', 'Format', '%s%s', 'ReadVariableNames', false);
color_names = strtrim(T{:,1});
hex_str     = char(strtrim(T{:,2}));
color_rgb   = [hex2dec(hex_str(:,2:3)), hex2dec(hex_str(:,4:5)), hex2dec(hex_str(:,6:7))] ./ 255;

rgb = repmat(default_rgb, numel(names), 1);
for i = 1:numel(names)
    ix = find(strcmpi(color_names, strtrim(string(names{i}))), 1);
    if ~isempty(ix)
        rgb(i, :) = color_rgb(ix, :);
    end
end
end
