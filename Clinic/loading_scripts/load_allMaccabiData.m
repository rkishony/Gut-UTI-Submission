% auto-script. Runs by `load_or_create`. Can also run manually from the command line.

%% Set out_filepath:
clear
folders = get_folders();
out_filepath = fullfile(folders.mat_outputs, 'allMaccabiData.mat');

%% Source Folders
fldr_translation_tables = [folders.translation_tables, filesep];
fldr_raw = [folders.Maccabi_data, filesep] ;
fldr_collection = [folders.collection, filesep] ;


%% get FOBT and UTI info from our collection
filenameFOBTtranslation='dam samui.xlsx';
filenameVTKtranslation='ids and req.xlsx';
filenamerFOBT='FOBTrecurrent.xlsx';
filename='FOBT_final.xlsx';

FOBT=readtable([fldr_collection filename],'Sheet','Our FOBT Collection','VariableNamingRule','preserve');
FOBT=extractCols(FOBT,[1,4,6],{'FOBT_ID','FOBT_Drisha','FOBT_Date'});
FOBT_Date=zeros(height(FOBT),1);
properdates=cellfun(@length,FOBT.FOBT_Date)==10;
FOBT_Date(properdates)=cellfun(@(x) datenum(x,'dd/mm/yyyy'), FOBT{properdates,'FOBT_Date'});
FOBT.FOBT_Date=FOBT_Date;

% Map FOBT_Drisha -> RandomID based on FOBT_dictionary
FOBT_dictionary=readtable([fldr_raw filenameFOBTtranslation]);
FOBT_dictionary=extractCols(FOBT_dictionary,[1,6],{'RandomID','FOBT_Drisha'});
FOBT = addMatch(FOBT, FOBT_dictionary, 'FOBT_Drisha' ,'RandomID');
FOBT(isnan(FOBT.RandomID),:) = [];
FOBT = sortrows(FOBT, {'RandomID', 'FOBT_Date'});
FOBT.FOBTisFirst = FOBT.RandomID ~= [-1; FOBT.RandomID(1:end-1)];

% Recurrent FOBT (not used)
rFOBT_dictionary=readtable([fldr_collection filenamerFOBT],'Sheet','Our collection');
rFOBT_dictionary=extractCols(rFOBT_dictionary,[1,2,4],{'FOBT_ID','FOBT_Drisha','RandomID'});
rFOBT_dictionary=addMatch(rFOBT_dictionary,FOBT,'RandomID','FOBT_ID','previousFOBT_ID');

UTIpop=readtable([fldr_collection filename],'Sheet','Scanning UTI','VariableNamingRule','preserve');
UTIpop=extractCols(UTIpop,[3,4],{'UTI_Drisha','UTI_ID'});

clear properdates  

%% Vitek dictionary to map Drisha -> RandomID
VTK_dictionary=readtable([fldr_raw filenameVTKtranslation]);
VTK_dictionary=extractCols(VTK_dictionary,[],{'RandomID','UTI_Drisha'});
VTK_dictionary=unique(VTK_dictionary,'rows');

% Deal with specific UTI_Drisha that are mapped to more than one RandomID:
fprintf('VTK_dictionary `UTI_Drisha` mapped to more than one `RandomID` (all  ): %d\n', height(VTK_dictionary)-height(unique(VTK_dictionary(:,'UTI_Drisha'))))
VTK_dictionary = VTK_dictionary(ismember(VTK_dictionary.RandomID, FOBT.RandomID), :);
fprintf('VTK_dictionary `UTI_Drisha` mapped to more than one `RandomID` (clean): %d\n', height(VTK_dictionary)-height(unique(VTK_dictionary(:,'UTI_Drisha'))))

%UTIpop = verbose_innerjoin(UTIpop, VTK_dictionary);
UTIpop = addMatch(UTIpop, VTK_dictionary, 'UTI_Drisha', 'RandomID');  % keeping unmatched lines

%% Load diagnosis
DIAG = readtable([fldr_raw 'diagnosis.xlsx']);
DIAG = extractCols(DIAG, [1,4,5], {'diagnosis_code','date_diagnosis','RandomID'});

% Identify UTI and bacteriuria
diagtranslation = load([fldr_translation_tables 'DiagnosisTableMarta.mat']).DiagnosisTable;
is_utis=contains(diagtranslation{:,'Description'},{'URINARY TRACT INFECTION','LOWER URINARY TRACT SYMPTOMS','LUTS','INFECTION URINARY','PYELONEPHRITIS','CYSTITIS'}) & strcmp(diagtranslation{:,'DiagnosticCategory'},'Genitourinary');
is_bu=contains(diagtranslation{:,'Description'},{'BACTERURIA'});
DIAG.UTI=ismember(DIAG.diagnosis_code,diagtranslation{is_utis,[1 2]});
DIAG.bacteruria=ismember(DIAG.diagnosis_code,diagtranslation{is_bu,[1 2]});

% Add specific description
[~,im1]=ismember(DIAG.diagnosis_code,diagtranslation.DiagnosisCode);
[~,im2]=ismember(DIAG.diagnosis_code,diagtranslation.ICD9);
im1(im1==0)=nan;
im2(im2==0)=nan;
im=min(im1,im2,'omitnan');
DIAG.description=repelem({'unknown'},height(DIAG))';
DIAG.description(im>0)=diagtranslation.Description(im(im>0));

% Convert date
DIAG.date_diagnosis=arrayfun(@(x) datenum(num2str(x),'yyyymmdd'),DIAG.date_diagnosis);

% Match with FOBT
DIAG = add_fobt(DIAG, FOBT, 'date_diagnosis');

%% Load Demographics
filename = 'demographics and rashmim.xlsx' ;
IDS=readtable([fldr_raw filename],'ReadVariableNames',false) ;
IDS = extractCols(IDS,[1,2,3,10,14,20],{'RandomID','DOB','Gender','Diabetes','Smoking','BMI'});
IDS.Properties.Description = filename ;
IDS.Gender = char(IDS.Gender) ; mf = ' MF' ; [~,im] = ismember(IDS.Gender,[1494,1504]) ; IDS.Gender = mf(im+1)' ;
IDS.age=2021-IDS.DOB;

% Match with FOBT:
IDS = add_fobt(IDS, FOBT, []);

%% Load Vitek (VTK, to map UTI_Drisha -> PhysicianID, SampleDate)
% TZNUM, TZ, LAB_RESP, DRISHA_NUM, PREG, DRISHA_COMPLETE, DATE_ISHUR, HAAVARA_DATE, PHSICIANID, TSHUVA_VITEK, HEARA_LINE, HEARA_TEXT_VITEK, SAMPLE_DATE, TEST_CODE
filename = 'vitek_titles.csv';
VTK = readtable([fldr_raw filename]) ;
VTK = extractCols(VTK, [4,7,8,9,13], {'UTI_Drisha','IshurDate','DrishaDate','PhysicianID','SampleDate'});
VTK = unique(VTK, 'rows');
VTK{:,{'IshurDate','DrishaDate','SampleDate'}} = yyyymmdd2datenum(VTK{:,{'IshurDate','DrishaDate','SampleDate'}});
assert(all(VTK{:,{'IshurDate'}} == VTK{:,{'DrishaDate'}}))
VTK(:, 'DrishaDate') = [];  % remove redundant dates
assert(all(VTK{:,{'SampleDate'}} <= VTK{:,{'IshurDate'}}))

%% Load Vitek mic (Bacteria identification and resistance)
filename='mic.csv' ;
warning off
mic = readtable([fldr_raw filename],'ReadVariableNames',true);
warning on
ttls=readtable([fldr_translation_tables 'MICfile_translation_all_2206.xlsx'],'Sheet','Field_names','VariableNamingRule','preserve') ;
ttls=ttls(ttls.In_Data>0,:);
mic=mic(:,ttls.field_name_in_data);
mic.Properties.VariableNames = ttls.Field_name;
%format them based on type
for i=1:height(ttls)
    fld = ttls.Field_name{i};
    switch ttls.Type{i}
        case {'YN','Letter'}, mic.(fld) = char(mic.(fld)) ;
        case 'Date', mic.(fld) = yyyymmdd2datenum(mic.(fld)) ;
    end
end

%% Create BAC
[BAC, ~, jBAC] = unique(mic(:,{'UTI_Drisha', 'Bacteria'})) ;
b = char(BAC.Bacteria) ; i = b(:,1)=='V' ; bn = zeros(size(i)) ; bn(i) = str2num(b(i,2:end)) ; BAC.Bacteria = bn ;

% Add resistances:
[RES_All,~,cant,kant] = sortunique(mic.AntCode) ;
RES_All.Properties.VariableNames = {'Code', 'Count'} ;

H = height(BAC);
L = height(RES_All);

cnts = accumarray([jBAC, kant], ones(numel(jBAC),1), [H, L]);
fprintf('Number of conlicting antibiotic calls: %g/%g (appearing in total of %g/%g bacteria)\n', ...
    sum(cnts(:)>1),sum(cnts(:)), sum(any(cnts>1,2)), sum(any(cnts,2)))

ind = sub2ind([H L], jBAC, kant);

v = char(zeros(H,L,'uint8')) ;
v(ind) = mic.Res;
BAC.v = v;

[~, vnum] = ismember(v, [char(0) 'SIR+']);
assert(all(vnum>0, 'all'))  % we only have S, I, R, + (and char(0) - no data)
vnum = vnum - 2;  % -1: NA,   0: S,   1: 'I',   2: 'R',   3: '+' (relevant only for ESBL)
vnum(vnum==-1) = nan;
BAC.vnum = vnum;

m = nan(H,L) ;
m(ind) = mic.MIC_val ;
BAC.mic = m;

sign = char(zeros(H,L,'uint8')) ;
sign(ind) = mic.MIC_sign ;
[~, sign] = ismember(sign, [char(0) '< >']);
assert(all(sign>0, 'all'))
sign = sign - 3;
sign(sign==-2) = nan;
BAC.sign = sign;

% Update BAC with PhysicianID, SampleDate (from VTK):
BAC = verbose_innerjoin(BAC, VTK);
BAC.SeniorHouse = floor(BAC.PhysicianID/100)==9999 & ~any(BAC.PhysicianID==[999985,999995,999980,999999,999996,999997],2) ; %%%
BAC(:, 'PhysicianID') = [];

% Use VTK_dictionary to get BAC's RandomID
BAC = verbose_innerjoin(BAC,VTK_dictionary);  % this reduce the number of lines

% Merge BAC with FOBTpop
BAC = add_fobt(BAC, FOBT, 'SampleDate');  % Same-patient FOBT lines are duplicated
BAC = sortrows(BAC,{'FOBT_ID','SampleDate'});

clear sign vnum v m H L mic

%% BAC -> SMP
% unique over all variables that are not isolate related (so we get sampes):
[SMP,iSMP,jSMP] = unique(BAC(:,{'UTI_Drisha', 'SampleDate', 'IshurDate', 'RandomID', 'FOBT_Date', 'FOBT_ID', 'FOBTisFirst', 'DateDiff'}), 'stable');
assert(height(SMP)==height(unique(BAC(:, {'UTI_Drisha', 'FOBT_ID'}))));

BAC.jSMP = jSMP;
[BAC.nBac, BAC.iBac] = j2ni(jSMP);
SMP.nBac = BAC.nBac(iSMP);

% Add SMP resistance (max over isolates).
vbac = BAC.vnum;
vsmp = nan(height(SMP), size(vbac,2));
for j = 1:height(SMP)
    rows = jSMP == j;
    if any(rows)
        vsmp(j,:) = max(vbac(rows,:),[],1);
    end
end
SMP.vnum = vsmp;

%% Add SampleDate from SMP to UTIpop
UTIpop = addMatch(UTIpop, SMP, 'UTI_Drisha', 'SampleDate');
UTIpop = add_fobt(UTIpop, FOBT, 'SampleDate', true);  % We keep missing lines

%% Load tested antibiotic names
filename = 'Antibiotics_full_new_paper_names_wGeneric.xls' ;
Abs = readtable([fldr_translation_tables filename]) ;
Abs.Properties.Description = filename ;
im = ismember(RES_All.Code,Abs.Code) ; 
codes = setdiff(RES_All.Code,Abs.Code)';
if ~isempty(codes)
    fprintf('\nCorrecting for antibiotic codes %s appearing in vitek file, but not in antibiotic list.\n', mat2str(codes'))
    warning off
    for code = codes
        Abs(end+1,{'Code', 'Name'}) = {code, sprintfc('Ab#%d',code)};
    end
    warning on
end
[~,j1] = ismember(RES_All.Code,Abs.Code) ;
[~,j2] = ismember(Abs.Code,RES_All.Code) ;
Abs = [Abs(j1,:); Abs(~j2,:)] ;
[~,j] = ismember(RES_All.Code, Abs.Code) ;
RES_All.Name = Abs.Name(j) ;

%% adding antibiotic purchases (not based on Largo);
filename = 'antibiotics.xlsx' ;
PRCH = readtable([fldr_raw filename],'ReadVariableNames',true)  ;
PRCH.Properties.Description = filename;
PRCH = extractCols(PRCH,[1 2 9 10 6 7],{'RandomID','Date','ATC_CD','ATC_Desc','Generic_CD','Generic_name'});
PRCH.Date = yyyymmdd2datenum(PRCH.Date) ;
PRCH = add_fobt(PRCH,FOBT);
PRCH = sortrows(PRCH, {'FOBT_ID', 'Date'});
PRCH.jPRCH = (1:height(PRCH))';

%% Load pathogen names
filename = 'Pathogens.xlsx' ;
Bugs = readtable([fldr_translation_tables filename]) ;
Bugs.Properties.Description = filename ;
Bugs.Code = cellfun(@(x) str2double(x(2:end)),Bugs.Code) ;
Bugs = concatinate_code_text(Bugs) ;
n = countoccurence(BAC.Bacteria,Bugs.Code) ; [n,k] = sort(n,'descend') ;
Bugs = Bugs(k,:) ;
Bugs = [Bugs; {0,'N/A'}] ;
codes = setdiff(BAC.Bacteria,Bugs.Code) ;
if ~isempty(codes)
    [tot,i] = countoccurence(BAC.Bacteria,codes) ;
    fprintf('\nSetting unknown bacterial codes %s to "Unknown".\n', mat2str(codes'))
    % disp(sortrows(table(cds,tot,'VariableNames',{'Code','Counts'}),-2))
end
Bugs = [Bugs; {-1,'Unknown'}] ;
[~,BAC.jBugs] = ismember(BAC.Bacteria,Bugs.Code) ;
BAC.jBugs(BAC.jBugs==0) = height(Bugs) ;

%% Antibiotic to ResFinder manual map
filename = 'antibiotics.txt' ;
Abs2ResFinder = readtable(fullfile(folders.meta_data, filename),'ReadVariableNames',true);
Abs2ResFinder = renamevars(Abs2ResFinder, 'Antibiotic', 'Name');
RES_All = addMatch(RES_All, Abs2ResFinder,'Name','ResFinder');

%% Save
save(out_filepath,'BAC','FOBT','UTIpop','DIAG','IDS','Abs','RES_All','SMP','PRCH','VTK','Bugs')

UTIclean = UTIpop(UTIpop.FOBTisFirst & ~isnan(UTIpop.FOBT_ID), {'UTI_ID', 'FOBT_ID', 'DateDiff'});
save(fullfile(folders.mat_outputs, 'UTIclean.mat'), 'UTIclean')


function tbl = add_fobt(tbl, fobt, date_fld, keep_missing_lines)
nA = inputname(1);
nB = inputname(2);
if nargin<3, date_fld='Date'; end
if nargin<4
    keep_missing_lines = false;
end
fobt = fobt(:, {'RandomID', 'FOBT_Date', 'FOBT_ID', 'FOBTisFirst'});
if keep_missing_lines
    tbl = addMatch(tbl, fobt, 'RandomID', {'FOBT_Date', 'FOBT_ID', 'FOBTisFirst'});
else
    tbl = verbose_innerjoin(tbl, fobt, 'names',{nA,nB});
end
if ~isempty(date_fld)
    tbl.DateDiff = tbl.(date_fld) - tbl.FOBT_Date;
end
end
