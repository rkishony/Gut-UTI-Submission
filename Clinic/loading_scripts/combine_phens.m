function [PHEN, PHEN_Entities] = combine_phens(PHEN1, PHEN2, PHEN1_Entities, PHEN2_Entities)
PHEN1 = PHEN1(PHEN1.rep==0, :);
PHEN2 = PHEN2(PHEN2.rep==0, :);

assert(isempty(setxor(PHEN2.FOBT_ID, PHEN1.FOBT_ID)));

PHEN1 = sortrows(PHEN1, 'FOBT_ID');
PHEN2 = sortrows(PHEN2, 'FOBT_ID');

assert(all(PHEN1.FOBT_ID == PHEN2.FOBT_ID))

PHEN = table;
PHEN.FOBT_ID = PHEN1.FOBT_ID;

PHEN_Entities = struct;

for fld = fieldnames(PHEN1_Entities)'
    f = fld{1};
    drugs_1 = cellfun(@(x)[x, '_A'], PHEN1_Entities.(f).Drug,'UniformOutput',false);
    drugs_2 = cellfun(@(x)[x, '_B'], PHEN2_Entities.(f).Drug,'UniformOutput',false);

    PHEN.(f) = [PHEN1.(f), PHEN2.(f)];
    PHEN_Entities.(f) = array2table([drugs_1; drugs_2], 'VariableNames', {'Drug'});
end
assert(all(PHEN1.isFOBT))
assert(all(PHEN2.isFOBT))
PHEN.isFOBT = true(height(PHEN), 1);

assert(all(~PHEN1.is_neg_ctrl))
assert(all(~PHEN2.is_neg_ctrl))
PHEN.is_neg_ctrl = false(height(PHEN), 1);

PHEN.isrep = false(height(PHEN), 1);

end
