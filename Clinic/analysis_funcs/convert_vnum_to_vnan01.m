function vnum = convert_vnum_to_vnan01(vnum, I_as_RS)
% Convert HMS vnum (nan, 0, 1, 2, 3) to vnan01 (nan 0 1). 
% Clasify 'I' as 'S'/'R'
%
% nan -> nan
% sen -> 0
% res -> 1
%

if I_as_RS == 'S'
    vnum_th = 2;  % I->S,  resistance is vnum>=2
else
    vnum_th = 1;  % I->R,  resistance is vnum>=1
end

vnum(~isnan(vnum)) = vnum(~isnan(vnum)) >= vnum_th;

end
