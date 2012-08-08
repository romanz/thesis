clc; 
clear;
D = dir('for_plot/*.mat');

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    S{k} = load(['for_plot/', f]);
end

S = reshape(S, 3, numel(S)/3);

betas = cellfun(@(s) s.betas(end), S);
[betas, I] = sort(betas, 2);
V = cellfun(@(s) s.V(end), S);

s = S{end, end};
sol = s.sol;
phi = regrid(sol.Phi);
r = sol.Phi.grid.r;
plot(convn(r, [1;1]/2, 'valid'), diff(phi) ./ repmat(diff(r), [1 size(phi, 2)]))