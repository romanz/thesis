function main_run(N, Rmax, betas)

%    N = 2.^(8);
%    Rmax = 100;
%    betas = 0.1 * 2.^[-2:0.5:5];
    for n=N 
        for r=Rmax
            g = grids(r, n+1, n+1, 1.0);
            sol = struct('alpha', 0.5, 'Du', 1.0, 'zeta', 10);
            fname = sprintf('sol_beta=%.3e_[%dx%d]_Rmax=%.1f_Du=%.2f_zeta=%.2f_alpha=%.2f.mat', ...
                betas(1), numel(g.r), numel(g.t), max(g.r), sol.Du, sol.zeta, sol.alpha);
            s = load(fname);
            init.Phi = regrid(s.sol.Phi);
            init.C   = regrid(s.sol.C);
            init.Vr  = regrid(s.sol.Vr);
            init.Vt  = regrid(s.sol.Vt);
            init.P   = regrid(s.sol.P);
            main_ephor(init, g, betas, sol);
        end
    end

end
