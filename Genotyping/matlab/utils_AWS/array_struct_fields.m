function A = array_struct_fields(S, flds)

for i = 1:numel(flds)
    fld = flds{i};
    if i == 1
        A = S.(fld);
    else
        A(i) = S.(fld);
    end
end
