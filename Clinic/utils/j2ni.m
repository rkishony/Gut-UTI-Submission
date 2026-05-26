function [n,i] = j2ni(x)
[xs,xk] = sort(x) ;
[~,xkk] = sort(xk) ;
[~,ui,uj] = unique(xs) ;
nt = hist(uj,1:length(ui))' ;
n = nt(uj) ;
i = (1:length(x))' ;
i = i-i(ui(uj))+1 ;
n = n(xkk) ;
i = i(xkk) ;
end
