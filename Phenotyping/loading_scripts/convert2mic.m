function mics = convert2mic(measurements, Abs, concs, nPlates)
Cons = concs.Concentration;

conv_vec = log2(Cons)-(nPlates+1);

mics = -1*ones(size(measurements));

for ab = 1:numel(Abs)
    conv = conv_vec(find(ismember(concs.Antibiotic, Abs{ab})));
    mics(:,ab) = 2.^(measurements(:,ab) + conv); 
end

end