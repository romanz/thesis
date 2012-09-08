function [init, g] = load_init(matfile)

S = load(matfile);
sol = S.sol;
g = S.g;
init.Phi = regrid(sol.Phi);
init.C   = regrid(sol.C);
init.Vr  = regrid(sol.Vr);
init.Vt  = regrid(sol.Vt);
init.P   = regrid(sol.P);
