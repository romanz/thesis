clear;
load 2012_08_11/100/256x256_Rmax=100_beta=0.100.mat
phi = regrid(sol.Phi);
r = sol.Phi.grid.r;
plot(convn(r, [1;1]/2, 'valid'), diff(phi) ./ repmat(diff(r), [1 size(phi, 2)]), '.-')
