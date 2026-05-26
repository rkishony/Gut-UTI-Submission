% Assert-only checks for deterministic SVG export environment.
% Run this at MATLAB startup/session begin:
%   run('from_dropbox2/assert_svg_export_env.m')

fprintf('Checking SVG export environment...\n');

% Expected values for this headless server setup.
expectedScreenPixelsPerInch = 72;
expectedScreenSize = [1 1 1024 768];
expectedMonitorPositions = [1 1 1024 768];
expectedDefaultFigurePosition = [232 246 560 420];
expectedDefaultAxesPosition = [0.13 0.11 0.775 0.815];
expectedDefaultAxesLooseInset = [0.13 0.11 0.095 0.075];
expectedDefaultAxesFontName = 'Helvetica';
expectedDefaultAxesFontSize = 10;
expectedDefaultTextFontName = 'Helvetica';
expectedDefaultTextFontSize = 10;
tol = 1e-10;

% DISPLAY must be set to a valid X11 display for correct Qt font metrics.
% If DISPLAY is stale (from old SSH session), restart MATLAB with current DISPLAY.
dispVal = getenv('DISPLAY');
assert(~isempty(dispVal), ...
    'DISPLAY is not set. MATLAB needs a valid X display for correct SVG font metrics.');
try
    % Verify the display is actually reachable
    f_test = figure('Visible','off');
    delete(f_test);
catch
    error('DISPLAY=%s is stale/invalid. Exit MATLAB and re-run matlab_tmux.', dispVal);
end

% Environment checks
assert(~usejava('desktop'), ...
    'Expected headless MATLAB: usejava(''desktop'') must be false.');

f_tmp = figure('Visible','off');
ax_tmp = axes(f_tmp);
rInfo = rendererinfo(ax_tmp);
delete(f_tmp);
fprintf('  Renderer: %s\n', rInfo.GraphicsRenderer);

assert(abs(get(groot, 'ScreenPixelsPerInch') - expectedScreenPixelsPerInch) <= tol, ...
    'ScreenPixelsPerInch mismatch.');

assert(all(abs(get(groot, 'ScreenSize') - expectedScreenSize) <= tol), ...
    'ScreenSize mismatch.');

assert(all(abs(get(groot, 'MonitorPositions') - expectedMonitorPositions) <= tol), ...
    'MonitorPositions mismatch.');

% Root defaults that influence figure geometry/text layout
assert(all(abs(get(groot, 'DefaultFigurePosition') - expectedDefaultFigurePosition) <= tol), ...
    'DefaultFigurePosition mismatch.');

assert(all(abs(get(groot, 'DefaultAxesPosition') - expectedDefaultAxesPosition) <= tol), ...
    'DefaultAxesPosition mismatch.');

assert(all(abs(get(groot, 'DefaultAxesLooseInset') - expectedDefaultAxesLooseInset) <= tol), ...
    'DefaultAxesLooseInset mismatch.');

assert(strcmp(get(groot, 'DefaultAxesFontName'), expectedDefaultAxesFontName), ...
    'DefaultAxesFontName mismatch.');

assert(abs(get(groot, 'DefaultAxesFontSize') - expectedDefaultAxesFontSize) <= tol, ...
    'DefaultAxesFontSize mismatch.');

assert(strcmp(get(groot, 'DefaultTextFontName'), expectedDefaultTextFontName), ...
    'DefaultTextFontName mismatch.');

assert(abs(get(groot, 'DefaultTextFontSize') - expectedDefaultTextFontSize) <= tol, ...
    'DefaultTextFontSize mismatch.');

% Font availability check (important for cross-viewer SVG text alignment)
availableFonts = listfonts;
assert(any(strcmpi(availableFonts, 'Helvetica')), ...
    'Helvetica not found in listfonts.');

fprintf('OK: SVG export environment matches expected settings.\n');
