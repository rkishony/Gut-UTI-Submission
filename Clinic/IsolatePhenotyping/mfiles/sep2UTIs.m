function BACres2 = sep2UTIs(MHSall)

%%%%% THIS SCRIPT acts on specific UTI IDs, ones that appear more than once
%%%%% in the maccabi data, belonging to the same patient, Maccabi's way of indicating a 
%%%%% mixed urine culture (multi-strain infection).
%%%%%
%%%%% This code is meant to change one of each pair of them into a different ID, the one
%%%%% that we gave it in the lab, in order to distinguish the two and avoid mixing them up.
%%%%% 
%%%%% The UTI rows are searched for, then I manually look at the maccabi-measured profile
%%%%% and compare it with what we measured in the lab. Finally, one is renamed according 
%%%%% to the profile matching.

ABs = {'Ampicillin'; 'Cefazolin'; 'Ceftazidime'; 'Ciprofloxacin'; 'Phosphomycin'; 'Nitrofurantoin'; 'Trimethoprim'};
Abs = {'AMP';'CEF';'CIP';'FOS';'NIT';'TMP'}; % so far based on as many ABs as possible.

UTIpop = MHSall.UTIpop;
cBACres = MHSall.BACres;

isAMP = cellfun(@(c) strcmp('Ampicillin',c), MHSall.RES.Name);
% isAMX = cellfun(@(c) strcmp('Amoxicillin-CA',c), MHSall.RES.Name);
isCEF = cellfun(@(c) strcmp('Cefazolin',c), MHSall.RES.Name); % I think if we want this to also work for cefalexin or otherwise
% isCEFT = cellfun(@(c) strcmp('Ceftazidime',c), MHSall.RES.Name); % I think if we want this to also work for cefalexin or otherwise
isCIP = cellfun(@(c) strcmp('Ciprofloxacin',c), MHSall.RES.Name); 
isFOS = cellfun(@(c) strcmp('Fosfomycin',c),MHSall.RES.Name); 
isNIT = cellfun(@(c) strcmp('Nitrofurantoin',c), MHSall.RES.Name); 
isTMP = cellfun(@(c) strcmp('Trimethoprim-Sulfa',c), MHSall.RES.Name); 
relABs = [find(isAMP) find(isCEF) find(isCIP) find(isFOS) find(isNIT) find(isTMP)]; 
% relABs=[14 16 3 2 13 7 6];

%% Make changes:

[G, UTI_Drisha_vals] = findgroups(cBACres.UTI_Drisha);
counts = splitapply(@(x) numel(unique(x)), cBACres.UTI_ID, G);
ambiguous = UTI_Drisha_vals(counts > 1);

if ~isempty(ambiguous)
    fprintf("Ambiguous UTI Drisha's detected, requiring manual resolution:\n")
    disp(cBACres(ismember(cBACres.UTI_Drisha, ambiguous),:))
end

match = [1 2; 1 2];
cBACres2 = cBACres(~ismember(cBACres.UTI_Drisha, ambiguous),:);
for i = 1:length(ambiguous)
    a = ambiguous(i);
    r = MHSall.BAC(MHSall.BAC.UTI_Drisha==a,:);
    % cur_rel_prof = r.vnum(:,relABs);
    % by eye, comparing cur_rel_prof = r.vnum(:,relABs);
    % with  IsolateResDataTable(ismember(IsolateResDataTable.uti_id,UTIpop.UTI_ID(UTIpop.UTI_Drisha==a)),:)
    t_to_add = r(match(i,:),:);
    t_to_add.UTI_ID = UTIpop.UTI_ID(UTIpop.UTI_Drisha==a);
    cBACres2 = [cBACres2; t_to_add];
    
end
% nonAmbiguousMask = ~ismember(UTIpop.UTI_Drisha, ambiguous);
% BACres2 = unique(BACres2(nonAmbiguousMask,:));

% the ones that need fixing are 466/469 and 623/627, but others which have
% the same problem are 257/258 belonging to fobt 257 and 630/636 belonging
% to fobt 550, which are not in the metagenomics dataset.

% relUTIs = [466; 257; 623; 630]; %; 580; 51]; % apparently the last two were our addition and we didn't get an extra sample from MHS
% replacewith = [469; 258; 627; 636]; %; 587; 52];
% ind_of_choice = [2; 1; 2; 2]; % was chosen at every iteration initially.
% % the last two were chosen mostly based on the bacterium species.
% for u = 1:length(relUTIs)
%     inds = find(ismember(cBACres.UTI_ID,relUTIs(u))); % there should be two
%     theprofiles = cBACres.vnum(inds,relABs)
% 
%     BACres2.UTI_ID(inds(ind_of_choice(u))) = replacewith(u);
% end
BACres2 = cBACres2;
end
