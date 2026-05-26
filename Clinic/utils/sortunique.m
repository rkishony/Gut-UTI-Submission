function [t,u,n,ind] = sortunique(s,N)
if nargin<2
    N = inf ;
end

if isnumeric(s)
    mynan = 12345 ;
    isn = isnan(s) ;
    s(isn) = mynan ; 
end

[u,~,j] = unique(s,'rows') ;
n = accumarray(j,ones(size(j)),[size(u,1),1]) ;
if isnumeric(s)
    u(u==mynan) = nan ;
end

[~,ns] = sort(n,'descend') ;
u = u(ns(1:min(N,end)),:) ;
n = n(ns(1:min(N,end))) ;
[~,ind] = ismember(s,u,'rows') ;
if istable(u)
    t = u ;
    t.Counts = n ;
else
    t = table(u,n,'VariableNames',{'Category','Counts'}) ;
end
end