function res = test_quadratic(init)
    % variables with boundary conditions
    % initial solution value
    grid = grids(logspace(0, 5, 300), linspace(0, pi, 80));
    if nargin < 1
        init.Phi = zeros(grid.Phi.size);
        init.C = ones(grid.C.size);
        init.Vr = zeros(grid.Vr.size);
        init.Vt = zeros(grid.Vt.size);
        init.P = zeros(grid.P.size);
    end
    sol = Solution(grid, init);
    
    sol.alpha = 0;
    sol.beta = 0.1;
    sol.gamma = 1.5;
    sol.Vinf = 0.1;
    
    [sol.bnd, sol.I] = boundary_conditions(sol);
    sol.eqn = system_equations(sol);
    
    v = sol.var;
    n = v.grid.numel;
    op = sol.eqn;
    r0 = op.res();
    G = op.grad();

    rng = RandStream('mt19937ar', 'Seed', 1);
    dv = 1e-5*rng.randn(n, 1);
    v.update(dv);
    
    r1 = op.res();
    dr = r1 - r0;
    ds = G*dv;
    plot([dr-ds])
    save main
    return;
end
