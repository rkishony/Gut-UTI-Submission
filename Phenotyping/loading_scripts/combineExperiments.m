function all_res = combineExperiments()

expnames = ["IsolateResA"    "IsolateResB"    "IsolateResC"];

data_fldr = fullfile(get_folders('source_data'), 'Phenotyping');
output_path = fullfile(get_folders('source_data'), 'Phenotyping', 'IsolateResData.csv');

distinguishing_conc = readtable(fullfile(get_folders('metadata'), 'distinguishing_conc.xlsx'));
Abs = distinguishing_conc.Properties.VariableNames;
AbsMIC = cellfun(@(x)[x '_MIC'], Abs, 'uni',0);
AllAbs_Growth = cellfun(@(x)['AllAbs_Growth_' x], Abs, 'uni',0);

NUM_FEATURES = 9;
REL_GROWTH_THRESHOLD = 0.2;
GROWTH_THRESHOLD = 0.02;

%% load resistance data

for exp = 1:numel(expnames)
    ename = expnames(exp);
    fname = fullfile(data_fldr, ename, 'results.csv');
    raw = readtable(fname);
    raw.expname = repmat(string(ename), height(raw), 1);
    raw = raw(:,[end, 1:end-1]);
    raw_measurements.(ename) = raw;
end

for exp = 1:2
    ename = expnames(exp);
    meas = raw_measurements.(ename);
    vars = meas.Properties.VariableNames(1:NUM_FEATURES+1);
    mics = meas{:, AbsMIC};
    res01 = double(mics > distinguishing_conc{:,:});
    res01(isnan(mics)) = nan;
    measurements.(ename) = [meas(:, vars), array2table(res01, 'VariableNames', Abs)];
end

for exp = 3
    ename = expnames(exp);
    meas = raw_measurements.(ename);
    vars = meas.Properties.VariableNames(1:NUM_FEATURES+1);
    growth = meas{:, AllAbs_Growth};
    growth0 = meas.AllAbs_Growth_MC;
    norm_growth = growth ./ growth0;
    res01 = double(growth >= GROWTH_THRESHOLD & norm_growth >= REL_GROWTH_THRESHOLD);
    res01(growth0 < GROWTH_THRESHOLD, :) = nan;
    measurements.(ename) = [meas(:, vars), array2table(res01, 'VariableNames', Abs)];
end

all_res = cat(1, measurements.IsolateResA, measurements.IsolateResB, measurements.IsolateResC);

writetable(all_res, output_path);

end