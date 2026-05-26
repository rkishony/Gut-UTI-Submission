function [PHEN, PHEN_Entities] = process_PHEN(PHEN, PHEN_Abs, max_col, min_col)
% Normalize the meta-resistance measurements by the MC control.

% Remove counts in `merged` plating sports:
counts = PHEN.counts;

% Assign 'lawns' -> max_col:
counts(isinf(PHEN.counts)) = max_col;
PHEN.clean_counts = counts;

% Remove samples with too few counts on the MC control:
counts(counts(:,1) < min_col,:) = nan;

% Normalize by MC:
PHEN.clean_counts_n = counts ./ counts(:,1);

% Define PHEN_Entities listing the drug names for each of the counts
% columns of PHEN:
drug_tbl = cell2table(PHEN_Abs', 'VariableNames', {'Drug'});
PHEN_Entities.counts = drug_tbl;
PHEN_Entities.clean_counts = drug_tbl;
PHEN_Entities.clean_counts_n = drug_tbl;
PHEN_Entities.ismerged = drug_tbl;

end