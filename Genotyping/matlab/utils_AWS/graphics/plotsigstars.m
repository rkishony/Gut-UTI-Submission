function [h_line, h_text] = plotsigstars(xspan, y, txt, vert)
% Draw a significance bracket between x0 and x1 above bars whose top is at
% y
% x0, x1 : xaapositions of the two bars
% y      : yaavalue of the bar tops
% txt    : label (e.g. '*', 'n.s.')

if nargin<4
    vert = 0.01;
end

hold on;
yrng = diff(ylim);
dh   = vert * yrng; % bracket height
y0   = y;
y1   = y + dh;

h_line = line(xspan([1 1 2 2]), [y0 y1 y1 y0], 'Color','k');
h_text = text(mean(xspan), y1 + 0.5*dh, txt, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom');
end
