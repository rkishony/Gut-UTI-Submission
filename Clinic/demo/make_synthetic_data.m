function make_synthetic_data(out_file)
% MAKE_SYNTHETIC_DATA  Generate a small, fully synthetic dataset for the demo.
%
% This creates ARTIFICIAL data with the same structure as the analysis-ready
% tables used by the Clinic pipeline (the real study data are protected
% patient records from Maccabi Health Services and cannot be shared).
%
% A genotype/phenotype -> urine-culture-resistance association is deliberately
% built in so that the demo figures are non-degenerate and illustrate the
% analysis. The numbers are random and carry no biological meaning.
%
% Output: a .mat file with the variables consumed by run_demo.m:
%   FOBT, BAC, SMP, IDS, DIAG, UTIpop, RES_All,
%   RESFDR, RESFDR_COV, RESFDR_MUT, PHEN, PHEN_Entities

if nargin < 1
    here = fileparts(mfilename('fullpath'));
    out_file = fullfile(here, 'demo_data.mat');
end

rng(1);  % reproducible

nPat  = 180;                 % number of patients (one first FOBT each)
day0  = datenum('2019-01-01');

% ---- Antibiotics / resistance reference table (RES_All) -----------------
% vnum columns correspond to these rows (>=20 rows required by the pipeline).
names = {'Ampicillin','Cefazolin','Ceftazidime','Cefuroxime axetil', ...
         'Ciprofloxacin','Fosfomycin','Nitrofurantoin','Trimethoprim_Sulfa', ...
         'Gentamicin','Sulfamethoxazole', ...
         'Amikacin','Meropenem','Piperacillin','Amoxicillin','Cephalexin', ...
         'Doxycycline','Levofloxacin','Norfloxacin','Tobramycin','Ertapenem'}';
resfinder = {'beta-lactam','beta-lactam','beta-lactam','beta-lactam', ...
             'quinolone','fosfomycin','nitrofurantoin','trimethoprim', ...
             'aminoglycoside','sulphonamide', ...
             'other','other','other','other','other', ...
             'other','other','other','other','other'}';
nAb = numel(names);
RES_All = table((1:nAb)', zeros(nAb,1), names, resfinder, ...
    'VariableNames', {'Code','Count','Name','ResFinder'});

% Column indices (into vnum) for the antibiotic classes we model:
col_beta = [1 2 3 4];   % beta-lactam / cephalosporins
col_cpr  = 5;           % ciprofloxacin (quinolone)
col_fos  = 6;           % fosfomycin
col_nit  = 7;           % nitrofurantoin
col_tmp  = 8;           % trimethoprim/sulfa

% ---- Per-patient latent "resistance" state per class --------------------
% G.<class>(p) = true means patient p carries resistance in that class.
classes = {'beta','cpr','fos','nit','tmp'};
G = struct();
for c = 1:numel(classes)
    G.(classes{c}) = rand(nPat,1) < 0.4;
end

% ---- Demographics (IDS) -------------------------------------------------
RandomID = (1:nPat)';
gender = repmat("M", nPat, 1);
gender(rand(nPat,1) < 0.6) = "F";
age = randi([25 92], nPat, 1);
% ensure a populated Female 50-80 subgroup:
force_f5080 = randperm(nPat, 45);
gender(force_f5080) = "F";
age(force_f5080) = randi([50 80], numel(force_f5080), 1);
IDS = table(RandomID, gender, age, 'VariableNames', {'RandomID','Gender','age'});

% ---- FOBT (faecal samples) ---------------------------------------------
FOBT_ID   = RandomID;                       % one first FOBT per patient
FOBT_Date = day0 + randi(700, nPat, 1);
FOBTisFirst = true(nPat,1);
FOBT = table(RandomID, FOBT_ID, FOBT_Date, FOBTisFirst);

% ---- Urine cultures (SMP + BAC) ----------------------------------------
% Each patient gets 1-4 positive urine cultures within [0,365] days.
smp_rows = {};   % accumulate structs
bac_rows = {};
uti_drisha = 0;
species_main = 512;   % pretend E. coli code
for p = 1:nPat
    nP = randi([1 4]);
    % first culture is "coinciding" (0-3 days) for ~half the patients
    dd = [];
    if rand < 0.5
        dd(end+1) = randi([0 3]);
    end
    while numel(dd) < nP
        dd(end+1) = randi([5 360]); %#ok<AGROW>
    end
    dd = sort(dd);
    for k = 1:numel(dd)
        uti_drisha = uti_drisha + 1;
        % resistance calls per antibiotic column (nan = not measured):
        v = nan(1, nAb);
        v(:) = nan;
        measured = rand(1, nAb) < 0.9;
        pres_prob = 0.2 * ones(1, nAb);
        if G.beta(p), pres_prob(col_beta) = 0.6; end
        if G.cpr(p),  pres_prob(col_cpr)  = 0.6; end
        if G.fos(p),  pres_prob(col_fos)  = 0.6; end
        if G.nit(p),  pres_prob(col_nit)  = 0.6; end
        if G.tmp(p),  pres_prob(col_tmp)  = 0.6; end
        modelled = false(1,nAb); modelled([col_beta col_cpr col_fos col_nit col_tmp]) = true;
        callcols = find(measured & modelled);
        v(callcols) = 2 * (rand(1, numel(callcols)) < pres_prob(callcols)); % 0=S, 2=R
        % occasional species difference from the coinciding culture:
        sp = species_main;
        if rand < 0.15, sp = species_main + randi(5); end

        s = struct('UTI_Drisha', uti_drisha, 'SampleDate', FOBT_Date(p)+dd(k), ...
            'RandomID', p, 'FOBT_Date', FOBT_Date(p), 'FOBT_ID', p, ...
            'FOBTisFirst', true, 'DateDiff', dd(k), 'nBac', 1);
        s.vnum = v;
        smp_rows{end+1} = s; %#ok<AGROW>

        b = s;
        b.Bacteria = sp;
        bac_rows{end+1} = b; %#ok<AGROW>
    end
end

SMP = struct2table_vnum(smp_rows, nAb);
BAC = struct2table_vnum(bac_rows, nAb);

% Sort SMP by FOBT_ID then DateDiff (required by merge/addClosest):
SMP = sortrows(SMP, {'FOBT_ID','DateDiff'});
BAC = sortrows(BAC, {'FOBT_ID','DateDiff'});

% ---- Diagnoses (DIAG) ---------------------------------------------------
% Give ~60% of urine cultures a coinciding UTI diagnosis.
di_fobt = []; di_date = []; di_uti = [];
for i = 1:height(SMP)
    if rand < 0.6
        di_fobt(end+1,1) = SMP.FOBT_ID(i);        %#ok<AGROW>
        di_date(end+1,1) = SMP.SampleDate(i) + randi([-2 2]); %#ok<AGROW>
        di_uti(end+1,1)  = true;                  %#ok<AGROW>
    end
end
DIAG = table(di_fobt, di_date, logical(di_uti), false(size(di_uti)), ...
    'VariableNames', {'FOBT_ID','date_diagnosis','UTI','bacteruria'});

% ---- UTIpop (isolate collection bookkeeping; minimal) ------------------
UTIpop = table((1:height(SMP))', SMP.FOBT_ID, true(height(SMP),1), SMP.DateDiff, ...
    'VariableNames', {'UTI_ID','FOBT_ID','FOBTisFirst','DateDiff'});

% ---- ResFinder sample table (RESFDR) -----------------------------------
% One metagenome per patient FOBT.
Sample  = arrayfun(@(x) sprintf('F%d',x), FOBT_ID, 'UniformOutput', false);
bpTotal = 1.5e9 * ones(nPat,1);      % avg genome coverage ~300x -> isGood
ReadsTotal = round(bpTotal/300);
RESFDR = table(Sample, true(nPat,1), FOBT_ID, false(nPat,1), ReadsTotal, bpTotal, ...
    (1:nPat)', 'VariableNames', ...
    {'Sample','isFOBT','FOBT_ID','isrep','ReadsTotal','bpTotal','jRF'});

% ---- ResFinder gene coverage (RESFDR_COV) ------------------------------
% Alleles grouped by Drug (resfinder class) and Gene.
% Depth is an absolute k-mer depth; the pipeline normalises by avg coverage
% (bpTotal/5e6 ~ 300) to get normalised gene coverage cov_Drug.
avgcov = bpTotal(1) / 5e6;
alleles = {  % Drug(class)      Gene       Allele
    'beta_lactam',  'blaTEM',  'blaTEM-1'
    'beta_lactam',  'blaCTX',  'blaCTX-M-15'
    'beta_lactam',  'blaSHV',  'blaSHV-12'
    'fosfomycin',   'fosA',    'fosA3'
    'quinolone',    'qnrS',    'qnrS1'
    'quinolone',    'qnrB',    'qnrB19'
    'trimethoprim', 'dfrA',    'dfrA17'
    'sulphonamide', 'sul',     'sul1'
    'aminoglycoside','aac',    'aac(3)-IIa'
    };
allele_class = { ... % which latent class controls each allele row
    'beta','beta','beta','fos','cpr','cpr','tmp','tmp','beta'};

cS = {}; cDrug = {}; cGene = {}; cAllele = {}; cDepth = [];
for p = 1:nPat
    for a = 1:size(alleles,1)
        cls = allele_class{a};
        isRes = G.(cls)(p);
        if isRes
            covval = 0.1 + 0.85*rand;      % resistant: high normalised coverage
        else
            if rand < 0.5, continue; end   % often simply absent
            covval = 0.001 + 0.03*rand;    % otherwise trace coverage
        end
        cS{end+1,1}     = RESFDR.Sample{p};   %#ok<AGROW>
        cDrug{end+1,1}  = alleles{a,1};       %#ok<AGROW>
        cGene{end+1,1}  = alleles{a,2};       %#ok<AGROW>
        cAllele{end+1,1}= alleles{a,3};       %#ok<AGROW>
        cDepth(end+1,1) = covval * avgcov;    %#ok<AGROW>
    end
end
nCov = numel(cS);
RESFDR_COV = table(cS, cDrug, cGene, cAllele, cDepth, ...
    99*ones(nCov,1), 99*ones(nCov,1), false(nCov,1), ...
    'VariableNames', {'Sample','Drug','Gene','Allele','Depth', ...
                      'Template_Coverage','Query_Identity','Singleton'});

% ---- ResFinder point mutations (RESFDR_MUT) ----------------------------
% Quinolone-resistance mutations only (single Drug value -> single mut_Drug col).
mS = {}; mDrug = {}; mGene = {}; mAllele = {};
mut_alleles = {'gyrA','gyrA p.S83L'; 'gyrA','gyrA p.D87N'; 'parC','parC p.S80I'};
for p = 1:nPat
    if G.cpr(p)
        nMut = randi([2 3]);
    elseif rand < 0.1
        nMut = 1;
    else
        nMut = 0;
    end
    for m = 1:nMut
        mS{end+1,1}     = RESFDR.Sample{p};        %#ok<AGROW>
        mDrug{end+1,1}  = 'Ciprofloxacin';         %#ok<AGROW>
        mGene{end+1,1}  = mut_alleles{m,1};        %#ok<AGROW>
        mAllele{end+1,1}= mut_alleles{m,2};        %#ok<AGROW>
    end
end
nMutTot = numel(mS);
RESFDR_MUT = table(mS, mDrug, mGene, mAllele, ones(nMutTot,1), false(nMutTot,1), ...
    'VariableNames', {'Sample','Drug','Gene','Allele','Depth','Singleton'});

% ---- Community phenotyping (PHEN + PHEN_Entities) ----------------------
% Combined table as produced by process_PHEN + combine_phens: normalised
% fraction of resistant colonies per antibiotic (MC control == 1).
phen_pat = sort(randperm(nPat, 140))';   % subset with phenotyping
drugsA = {'MC_A','CEF_A','TMP_A','CPR_A'};      % experiment A (+_A suffix)
drugsB = {'MC_B','NIT-8_B','FOS-16_B'};         % experiment B (+_B suffix)
phen_drugs = [drugsA, drugsB];
% class controlling each phen column (MC has none):
phen_class = {'', 'beta', 'tmp', 'cpr', '', 'nit', 'fos'};

nP = numel(phen_pat);
ccn = ones(nP, numel(phen_drugs));
for i = 1:nP
    p = phen_pat(i);
    for d = 1:numel(phen_drugs)
        cls = phen_class{d};
        if isempty(cls), ccn(i,d) = 1; continue; end   % MC control
        if G.(cls)(p)
            ccn(i,d) = 0.1 + 0.85*rand;    % resistant: high colony fraction
        else
            ccn(i,d) = 0.005 + 0.035*rand; % sensitive: low fraction
        end
    end
end
PHEN = table(FOBT_ID(phen_pat), 'VariableNames', {'FOBT_ID'});
PHEN.clean_counts_n = ccn;
PHEN.isFOBT     = true(nP,1);
PHEN.is_neg_ctrl= false(nP,1);
PHEN.isrep      = false(nP,1);

drug_tbl = cell2table(phen_drugs', 'VariableNames', {'Drug'});
PHEN_Entities = struct();
PHEN_Entities.counts        = drug_tbl;
PHEN_Entities.clean_counts  = drug_tbl;
PHEN_Entities.clean_counts_n= drug_tbl;
PHEN_Entities.ismerged      = drug_tbl;

% ---- Save ---------------------------------------------------------------
save(out_file, 'FOBT','BAC','SMP','IDS','DIAG','UTIpop','RES_All', ...
    'RESFDR','RESFDR_COV','RESFDR_MUT','PHEN','PHEN_Entities');
fprintf('Wrote synthetic demo data: %s\n', out_file);
fprintf('  Patients: %d | urine cultures: %d | phenotyped: %d\n', ...
    nPat, height(SMP), nP);

end


function T = struct2table_vnum(rows, nAb)
% Assemble a table from a cell array of scalar structs whose `vnum` field is
% a 1xnAb row (kept as a matrix column in the table).
n = numel(rows);
UTI_Drisha  = zeros(n,1);
SampleDate  = zeros(n,1);
RandomID    = zeros(n,1);
FOBT_Date   = zeros(n,1);
FOBT_ID     = zeros(n,1);
FOBTisFirst = false(n,1);
DateDiff    = zeros(n,1);
nBac        = zeros(n,1);
vnum        = nan(n, nAb);
hasBac      = isfield(rows{1}, 'Bacteria');
Bacteria    = zeros(n,1);
for i = 1:n
    r = rows{i};
    UTI_Drisha(i)  = r.UTI_Drisha;
    SampleDate(i)  = r.SampleDate;
    RandomID(i)    = r.RandomID;
    FOBT_Date(i)   = r.FOBT_Date;
    FOBT_ID(i)     = r.FOBT_ID;
    FOBTisFirst(i) = r.FOBTisFirst;
    DateDiff(i)    = r.DateDiff;
    nBac(i)        = r.nBac;
    vnum(i,:)      = r.vnum;
    if hasBac, Bacteria(i) = r.Bacteria; end
end
T = table(UTI_Drisha, SampleDate, RandomID, FOBT_Date, FOBT_ID, ...
    FOBTisFirst, DateDiff, nBac);
T.vnum = vnum;
if hasBac
    T.Bacteria = Bacteria;
end
end
