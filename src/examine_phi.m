clear;
load for_irad/alpha=0.5_Du=1_zeta=10_large_betas/100/256x256_Rmax=100_beta=5.800_alpha=0.5_Du=1.0_zeta=10.0.mat
phi = regrid(sol.Phi);
r = sol.Phi.grid.r;
plot(convn(r, [1;1]/2, 'valid'), diff(phi) ./ repmat(diff(r), [1 size(phi, 2)]), '.-')
