function zeta = zeta_potential(sol)
g = Grid(1, sol.grid.Vt.t(2:end-1));
zeta = sol.zeta - log(Interp(g, sol.C));
zeta = zeta.res;
