function fix_svg_size(infile, targetW, targetH, outfile, units)

if nargin<4, outfile = infile; end
if nargin<5, units = 'in'; end


%--- Read entire SVG into a string
txt = fileread(infile);

%--- Extract original numeric width & height
tok = regexp(txt, ...
        '<svg[^>]*width=["'']([\d.]+)["''][^>]*height=["'']([\d.]+)["'']', ...
            'tokens','once');
if isempty(tok)
        error('Couldn''t find width/height on the <svg> tag.');
end
origW = str2double(tok{1});
origH = str2double(tok{2});

%--- Build new attribute strings
newW = sprintf('%g%s', targetW, units);
newH = sprintf('%g%s', targetH, units);
viewBoxStr = sprintf('viewBox="0 0 %g %g"', origW, origH);

%--- Replace width in <svg>
txt = regexprep(txt, ...
    '(<svg[^>]*?)\swidth=["''][^"'']*["'']', ...
    ['$1 width="' newW '"'], ...
    'once');
%--- Replace height
txt = regexprep(txt, ...
    '(<svg[^>]*?)\sheight=["''][^"'']*["'']', ...
    ['$1 height="' newH '"'], ...
    'once');
%--- Replace or insert viewBox
if contains(txt, 'viewBox=')
    txt = regexprep(txt, ...
        '(<svg[^>]*?)\sviewBox=["''][^"'']*["'']', ...
        ['$1 ' viewBoxStr], ...
        'once');
else
    txt = regexprep(txt, ...
        '<svg', ...
        ['<svg ' viewBoxStr], ...
        'once');
end

%--- Write out
fid = fopen(outfile,'w');
if fid<0, error('Cannot open %s for writing.', outfile); end
fwrite(fid, txt, 'char');
fclose(fid);
end
