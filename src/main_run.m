function main_run(N, Rmax, betas)

%    N = 2.^(8);
%    Rmax = 100;
%    betas = 0.1 * 2.^[-2:0.5:5];
    for n=N 
        for r=Rmax
            g = grids(r, n+1, n+1, 1.0);
            main_ephor([], g, betas, struct('alpha', 0.5, 'Du', 1.0, 'zeta', 10));
        end
    end

end
