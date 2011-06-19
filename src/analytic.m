function sol = analytic(sol)
    grid = sol.grid;
    sol.Vinf = -sol.beta * 2 * log((1 + 1/sqrt(sol.gamma)) / 2);
    
    sol.Phi =  sol.beta  * (0.25 * grid.Phi.X.^(-2) - grid.Phi.X) .* cos(grid.Phi.Y);
    sol.C = 1 + sol.beta  * 0.75 * grid.C.X.^(-2)                 .* cos(grid.C.Y);
    sol.Vx =  sol.Vinf     * (1 - grid.Vx.X.^(-3))                .* cos(grid.Vx.Y);
    sol.Vy = -sol.Vinf * (1 + 0.5*grid.Vy.X.^(-3))                .* sin(grid.Vy.Y);
    sol.P = zeros(grid.P.sz);
end