clear
load kuk


%% Find matching PRCH for each SMP (if single match)

SMP.jPRCH = zeros(height(SMP), 1);
SMP.isEmpirical = false(height(SMP),1);
SMP.isSameDay = false(height(SMP),1);
for i = 1:height(SMP)
    isPRCHsameid = find(PRCH.FOBT_ID == SMP.FOBT_ID(i));
    assert(all(diff(PRCH.Date(isPRCHsameid)) >= 0))
    is_matched_date = inrange(PRCH.Date(isPRCHsameid) - SMP.SampleDate(i), [-1, 3]);
    jPRCH_sameid_and_matched_date = isPRCHsameid(is_matched_date);
    if numel(jPRCH_sameid_and_matched_date) == 0
        continue
    end
    if numel(jPRCH_sameid_and_matched_date) == 1 || diff(PRCH.Date(jPRCH_sameid_and_matched_date(1:2))) > 0
        chosen_jPRCH = jPRCH_sameid_and_matched_date(1);
        SMP.jPRCH(i) = chosen_jPRCH;
        SMP.isEmpirical(i) = PRCH.Date(chosen_jPRCH) < SMP.IshurDate(i);
        SMP.isSameDay(i) = PRCH.Date(chosen_jPRCH) <= SMP.SampleDate(i);
    end

end

SMP = addMatch(SMP, PRCH, 'jPRCH', 'ATC_Desc');
SMP = addMatch(SMP, PRCH, 'jPRCH', 'Date', 'DatePRCH');
SMP.PRCH_diff = SMP.DatePRCH - SMP.SampleDate;

%%
figure(1);clf
treatment = SMP.jPRCH > 0;
empirical = SMP.isEmpirical;
dt = SMP.DatePRCH - SMP.SampleDate;
dt_u = -3:4;
c1 = histc(dt(treatment & empirical), dt_u);
c2 = histc(dt(treatment & ~empirical), dt_u);
bar(dt_u, [c1 c2])
xlim([-1.5 3.5])
xlabel('Time difference between antibiotic purchase and UTI sample date')
legend({'Empirical', 'Not empirical'}, Box='off')

%% Define PRCH-to-RES_All match table

Match = cell2table({
    "amoxicillin"       "Amoxicillin-CA"
    "cefalexin"         "Cephalexin"  % "Cefazolin", "Ceftazidime"
    "cefuroxime"        "Cefuroxime axetil"
    "ciprofloxacin"     "Ciprofloxacin"
    "fosfomycin"        "Fosfomycin"
    "methenamine"       ""
    "nitrofurantoin"    "Nitrofurantoin"
    "ofloxacin"         "Ofloxacin"
    "rifaximin"         ""
    }, "VariableNames", ["ATC_Desc", "RES_All"]);

Match.jRES_All = cellfun(@(s) selectOutput(@ismember, 2, s, string(RES_All.Name)), Match.RES_All);
jRES = Match.jRES_All(Match.jRES_All>0);
[~, Match.jjRES] = ismember(Match.jRES_All, jRES);
Match.jMatch = (1:height(Match))';

%% RF Predictions
RF_drugs = ["beta_lactam", "fosfomycin", "TMP_SMX", "CPR"];
jBetaLactam = find(strcmp(RF_Entities.cov_Drug.Drug, "beta_lactam"));
jFOS = find(strcmp(RF_Entities.cov_Drug.Drug, "fosfomycin"));
RF.values = [RF.cov_Drug(:,[jBetaLactam, jFOS]), RF.cov_TMP_SMX, RF.mut_Drug];
SMP = addMatch(SMP, RF, 'FOBT_ID', 'values', 'RF_values');


%% PHEN Predictions
PHEN_drugs = ["CEF", "NIT", "TMP", "CIP"];
[~, j] = ismember(PHEN_drugs, PHEN_Entities.clean_counts_n.Drug);
PHEN.values = PHEN.clean_counts_n(:, j);
SMP = addMatch(SMP, PHEN, 'FOBT_ID', 'values', 'PHEN_values');


%%
SMP = addMatch(SMP, Match, 'ATC_Desc', {'jMatch', 'jjRES'});
vnum = SMP.vnum;
vnum(isnan(vnum)) = -1;


%% Get the resistance of the treated drug
treated = SMP.jjRES > 0;

SMP.PRCH_vnum = nan(height(SMP),1);
for i = find(treated)'
    SMP.PRCH_vnum(i) = vnum(i, jRES(SMP.jjRES(i)));
end


%% Numbers
measured = SMP.PRCH_vnum >= 0 & SMP.FOBTisFirst;
withinYear = inrange(SMP.DateDiff, [-1 365]);

okPHEN = PHEN.clean_counts(:,1) > 20;

[~, SMP.jRF] = ismember(SMP.FOBT_ID, RF.FOBT_ID);
[~, SMP.jokPHEN] = ismember(SMP.FOBT_ID, PHEN.FOBT_ID(okPHEN));

hasFOBTmeasurement = SMP.jRF>0 | SMP.jokPHEN>0;

fprintf('Total PUC treated and measured %d\n', sum(measured))
fprintf('Total PUC treated and measured within a year %d\n', sum(measured & withinYear))
fprintf('Total PUC treated and measured within a year and empirical %d\n', sum(measured & withinYear & SMP.isEmpirical))
fprintf('Total PUC treated and measured within a year and empirical and RF %d\n', sum(measured & withinYear & SMP.isEmpirical & hasFOBTmeasurement))
ok = measured & withinYear & SMP.isEmpirical & hasFOBTmeasurement;
wrong = SMP.PRCH_vnum > 0;
fprintf('Total Patients treated and measured and RF within a year %d\n', numel(unique(SMP.FOBT_ID(ok))))

fprintf('Fraction mistakes in empirical %d, %4.2f\n', sum(SMP.PRCH_vnum(ok)>0), mean(SMP.PRCH_vnum(ok)>0))






%

figure(2);clf
imagesc(vnum(measured, jRES))
hold on
plot(SMP.jjRES(measured), 1:sum(measured), 'xk');
