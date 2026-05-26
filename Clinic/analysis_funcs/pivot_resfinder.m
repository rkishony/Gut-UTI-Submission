function [RF, RF_Entities, MUTCOV] = pivot_resfinder(RF, MUTCOV, func, prefix, norm)
% Pivot over the coverage/mutation tables to create matrices of
% cov/mut per sample, per 'entity' (allele, gene, and drug).
%
% MUTCOV can be RESFDR_MUT or RESFDR_COV

if nargin<5
    norm = [];
end

[~, MUTCOV.jRF] = ismember(MUTCOV.Sample, RF.Sample);
MUTCOV(~MUTCOV.jRF,:) = [];

RF_EntNames = {'Drug', 'Gene', 'Allele'};
RF_Entities = struct;
for iEnt = 1:numel(RF_EntNames)
    genericEnt = RF_EntNames{iEnt};

    % pivot Samble X Entity
    Ent = [prefix genericEnt];
    RF_Entities.(Ent) = unique(MUTCOV(:, RF_EntNames(1:iEnt)), 'rows');
    tot = table_pivot(func, 0, MUTCOV, 'Depth', RF(:,'jRF'), RF_Entities.(Ent));
    if ~isempty(norm)
        tot = tot ./ norm;
    end
    RF.(Ent) = tot;
end

end
