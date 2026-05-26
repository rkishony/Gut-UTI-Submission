function res = get_resistance_matrix(T)
meta_vars = {'expname','plate_num','well','name','FOBT_ID','UTI_ID','isolate_id','rep','is_neg_ctrl','from_lawn'};
ab_vars = setdiff(T.Properties.VariableNames, meta_vars, 'stable');
res = T{:, ab_vars};
end
