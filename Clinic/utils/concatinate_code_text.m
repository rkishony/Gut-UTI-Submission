function ntbl = concatinate_code_text(tbl)
% code, num, text
tbl = sortrows(unique(tbl),[1 2]) ;
i1 = find(tbl{:,2}==1)' ;
ntbl = tbl(i1,[1 3]) ;
i1 = [i1,height(tbl)+1] ;
for i = 1:length(i1)-1
    txt = tbl{i1(i),3}{1} ;
    for j=i1(i)+1:i1(i+1)-1
        txt = [txt ' ' tbl{j,3}{1}] ;
    end
    ntbl{i,2} = {txt} ;
end
