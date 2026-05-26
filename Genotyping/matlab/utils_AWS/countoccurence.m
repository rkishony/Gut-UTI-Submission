function [n,im] = countoccurence(a,b) 
[ub,~,jb] = unique(b) ;
[im0,im] = ismember(a,ub) ;
nub = accumarray(im(im0),ones(sum(im0),1),[length(b),1]) ;
n = nub(jb) ;
end