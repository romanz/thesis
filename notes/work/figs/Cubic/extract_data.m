clc;
D = dir('*.mat')
for k = 1:numel(D)
    s = load(D(k).name, 'solutions');
    v{k} = cellfun(@(r) r.Vinf, s.solutions);
    b{k} = cellfun(@(r) r.beta, s.solutions);
    g{k} = cellfun(@(r) r.gamma, s.solutions);
end
save data.mat v b g
