function T = verbose_innerjoin(A, B, varargin)

if length(varargin)>=1 && strcmpi(varargin{1},'names')
    names = varargin{2};
    varargin = varargin(3:end);
else
    names = [];
end

kIdx = find(strcmpi(varargin,'Keys'),1);
if isempty(kIdx)
    joinKeys = intersect(A.Properties.VariableNames, B.Properties.VariableNames);
    keySource = 'auto';
    varargin = [varargin, {'Keys', joinKeys}];
else
    joinKeys = varargin{kIdx+1};
    if ischar(joinKeys)||isstring(joinKeys)
        joinKeys = cellstr(joinKeys);
    end
    keySource = 'explicit';
end

ka = groupcounts(A, joinKeys);      % has vars: joinKeys + GroupCount
kb = groupcounts(B, joinKeys);
ka = renamevars(ka, 'GroupCount', 'CountA');
kb = renamevars(kb, 'GroupCount', 'CountB');
dup = innerjoin(ka, kb, 'Keys', joinKeys);  % join on same keys
if any(dup.CountA > 1 & dup.CountB > 1)
    prod_msg = ', Cartesian products created!';
else
    prod_msg = '';
end

T = innerjoin(A, B, varargin{:});

if isempty(names)
    nA = inputname(1); if isempty(nA), nA = 'A'; end
    nB = inputname(2); if isempty(nB), nB = 'B'; end
else
    nA = names{1}; 
    nB = names{2}; 
end    
kStr = strjoin(joinKeys, ', ');
fprintf('[innerjoin] %s (%d rows) + %s (%d rows) -> %d rows on key(s): %s [%s]%s\n', ...
    nA, height(A), nB, height(B), height(T), kStr, keySource, prod_msg);
end
