function fobt_val = get_FOBT_values(DATA, Entities, entity, field)
% Get specified values from FOBT

if nargin<4
    field = "";
end

if field == ""
    j_field = 1;
    assert(size(DATA.(entity),2)==1)
else
    entities = Entities.(entity);
    j_field = find(strcmp(field, entities{:,end}));
end
fobt_val = DATA.(entity)(:,j_field);
if size(fobt_val,2) > 1
    fprintf('WARNING: Colapsing more than one fobt entities\n')
    fobt_val = max(fobt_val,[],2);
end
