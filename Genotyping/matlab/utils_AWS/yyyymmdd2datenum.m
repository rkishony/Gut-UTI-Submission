function d = yyyymmdd2datenum(x)

d = datenum([fix(x(:)/1e4),floor(mod(x(:),1e4)/1e2),mod(x(:),1e2)]) ;
d = reshape(d,size(x)) ;

return
