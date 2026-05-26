function [conf_mats, chosen, Fisher_p_values, KS_p_values] = ...
    get_fobt_puc_confusion_mat_and_association_statistics(...
    fobt_val, fobt_val_th, num_puc_measured, num_puc_res, nres)

% Select patients with sen only PUC, or more than 1:n_res_limit res:
chosen = [num_puc_res==0, num_puc_res>=(1:nres)];

% Restrict to having atleast one resistance measurement:
chosen = chosen & (num_puc_measured>=1);

% Remove nan fobt_val
chosen = chosen & ~isnan(fobt_val);

% Res vs Sen statistics:
fobt_val_sen = fobt_val(chosen(:,1));

conf_mats = nan(2 ,2, nres);
KS_p_values = nan(1, nres);
Fisher_p_values = nan(1, nres);

for i = 1:nres
    fobt_val_res = fobt_val(chosen(:,i+1));

    conf_mat = [
        sum(fobt_val_sen<=fobt_val_th), sum(fobt_val_res<=fobt_val_th)
        sum(fobt_val_sen>fobt_val_th), sum(fobt_val_res>fobt_val_th)
        ];
    conf_mats(:,:,i) = conf_mat;
    try
        [~, KS_p_values(i)] = kstest2(fobt_val_sen, fobt_val_res);
    catch
        KS_p_values(i) = nan;
    end
    [~, Fisher_p_values(i)] = fishertest(conf_mat); %%% add 'tail','right'
end
