function s12 = joinStruct(s1, s2)
s12 = s1; 
flds = fieldnames(s2);
for i = 1:numel(flds)
    fld = flds{i};
    s12.(fld) = s2.(fld);
end

end