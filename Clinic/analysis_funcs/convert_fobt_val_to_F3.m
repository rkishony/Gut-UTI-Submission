function F3 = convert_fobt_val_to_F3(fobt_val, fobt_th)
% Convert fobt value to 
% 1 nan
% 2 Sen
% 3 Res

F3 = 1 + (~isnan(fobt_val)) + (fobt_val > fobt_th);

end
