function [] = main
    % variables with boundary conditions
    [grid] = grids(logspace(0, 3, 50), linspace(0, pi, 20));
    
    % interior grid (ghost points removed)
    init.Phi = zeros(grid.Phi.size);
    init.C = ones(grid.C.size);
    init.Vr = zeros(grid.Vr.size);
    init.Vt = zeros(grid.Vt.size);
    init.P = zeros(grid.P.size);
    
    sol = Solution(grid, init);
    sol.alpha = 0;
    sol.beta = 0.1;
    sol.gamma = 1.5;
    sol.Vinf = 0.1;
    
    iters = 6;
    [bnd, I] = boundary_conditions(sol);
    P = select(~I)';
    Q = select(I)';
    [eqn] = equations(sol);
    [~, J] = find(sol.P.L); % Pressure indices
    for k = 1:iters
        for i = 1:5
            r = bnd.res();
            norm(r)
            if norm(r, inf) < 1e-14, break; end
            G = bnd.grad();
            H = G*P;
            dx = linsolve(H, r);
            sol.x.value = sol.x.value - P*dx;
        end
        G = bnd.grad();
        H = G*P; % select boundary
        G = linsolve(H, G);
        H = G(:, I);
        
        r = eqn.res();
        norm(r, inf)
        G = eqn.grad();
        G1 = G(:, I);
        G2 = G(:, ~I);
        G = G1 - G2 * H;
        dx = linsolve(G, r);
        sol.x.value = sol.x.value - Q*dx;
        sol.x.value(J) = sol.x.value(J) - mean(sol.x.value(J));
    end
end

function x = linsolve(A, B)
    sz = size(A);
    assert(sz(1) == sz(2));
    assert(sz(1) == size(B, 1))
    x = A \ B;
end

function [g] = grids(r, t)
    r = r(:);
    t = t(:);
    rg = [2*r(1) - r(2); r; 2*r(end) - r(end-1)];
    tg = [2*t(1) - t(2); t; 2*t(end) - t(end-1)];
    rc = average(rg, [1; 1]/2);
    tc = average(tg, [1; 1]/2);
    
    g.Phi = Grid(rc, tc);
    g.C = Grid(rc, tc);
    g.Vr = Grid(r, tc);
    g.Vt = Grid(rc, t);
    g.P = Grid(rc(2:end-1), tc(2:end-1));
    g.r = r;
    g.t = t;
end

function [bnd] = Boundary(op, dim, n)
    g = op.grid;
    switch dim
        case  1,    g = Grid(g.r(end-n+1:end), g.t(2:end-1));
        case -1,    g = Grid(g.r(1:n), g.t(2:end-1));
        case  2,    g = Grid(g.r(:), g.t(end-n+1:end));
        case -2,    g = Grid(g.r(:), g.t(1:n));
    end
    bnd = Selector(g, op); % get boundary
end

function [op, I] = boundary_conditions(sol)

    phi1 = Boundary(sol.Phi, -1, 2);
    logc1 = log(Boundary(sol.C, -1, 2));    
    
    g1 = Grid(sol.Vr.grid.r(1), sol.Vr.grid.t(2:end-1));
    bnd{1} = Interp(g1, phi1 + logc1); % Phi, C @ R=1
    bnd{2} = Deriv(g1, phi1 - logc1, 1); % Phi, C @ R=1
    
    g2 = Grid(sol.Vr.grid.r(end), sol.Vr.grid.t(2:end-1));
    phi2 = Boundary(sol.Phi, 1, 2);    
    field = @(r, t) -sol.beta*cos(t);
    bnd{3} = Deriv(g2, phi2, 1) - field; % Phi, C @ R=inf
    bnd{4} = Boundary(sol.C, 1, 1) - 1; % Phi, C @ R=inf
    
    bnd{5} = Boundary(sol.Vr, 1, 1) - ( @(r,t) -sol.Vinf*cos(t) );
    bnd{6} = Boundary(sol.Vt, 1, 1) - ( @(r,t)  sol.Vinf*sin(t) );
    
    bnd{7} = Symm(sol.Phi, 2);
    bnd{8} = Symm(sol.C, 2);
    bnd{9} = Symm(sol.Vr, 2);
    bnd{10} = Symm(sol.Vt, 1);
        
    g = Grid(1, sol.Vt.grid.t(2:end-1));
    zeta = -Interp(g, sol.Phi) - log(sol.gamma);
    
    phi   = Interp(Grid(1, sol.Phi.grid.t(2:end-1)), sol.Phi); % Phi at R=1
    Dphi  = Deriv(g, phi, 2); % Dphi/Dtheta at R=1 on Vt grid.
    Vs    = Func(zeta, '4*log((exp(x/2) + 1)/2)') * Dphi;
    
    bnd{11} = Interp(g, sol.Vt) - Vs; % Vt @ R=1
    bnd{12} = Boundary(sol.Vr, -1, 1); % Vr @ R=1
    
    op = Join(bnd{:});
    I = [interior(sol.Phi, 1); interior(sol.C, 1); ...
         interior(sol.Vr, 1); interior(sol.Vt, 1); interior(sol.P, 0)];
end

function I = interior(op, n)
    g = op.grid;
    I = false(g.size);
    I(1+n:end-n, 1+n:end-n) = true;
    I = I(:);
end

function eq = Symm(op, n)
    dim = 2;
    
    bnd = Boundary(op, +dim, n);
    if n == 1
        eq1 = bnd;
    else
        g = bnd.grid; g = Grid(g.r, mean(g.t));
        eq1 = Deriv(g, bnd, dim);
    end

    bnd = Boundary(op, -dim, n);    
    if n == 1
        eq2 = bnd;
    else
        g = bnd.grid; g = Grid(g.r, mean(g.t));
        eq2 = Deriv(g, bnd, dim);
    end
    
    eq = Join(eq1, eq2);
end

function flux = charge(sol)
    DPhi_Dr = Crop(Deriv(sol.Vr.grid, sol.Phi, 1), [0 1]);
    DPhi_Dt = Crop(Deriv(sol.Vt.grid, sol.Phi, 2), [1 0]);
    Cr = Interp(DPhi_Dr.grid, sol.C);
    Ct = Interp(DPhi_Dt.grid, sol.C);
    g = sol.Phi.grid;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2' * Cr * DPhi_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Ct * DPhi_Dt, 2) * '1/(r^2 * sin(t))';
    flux = fluxR + fluxT;
end

function flux = salt(sol)
    DC_Dr = Crop(Deriv(sol.Vr.grid, sol.C, 1), [0 1]);
    DC_Dt = Crop(Deriv(sol.Vt.grid, sol.C, 2), [1 0]);
    % XXX Add advection!
    g = sol.C.grid;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2' * DC_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * DC_Dt, 2) * '1/(r^2 * sin(t))';
    flux = fluxR + fluxT;
end

function flux = mass(sol)
    Dm_Dr = Crop(sol.Vr, [0 1]);
    Dm_Dt = Crop(sol.Vt, [1 0]);
    g = sol.P.grid;
    fluxR = Deriv(g, 'r^2' * Dm_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Dm_Dt, 2) * '1/(r * sin(t))';
    flux = fluxR + fluxT;
end

function [forceR, forceT] = momentum(sol)
    Er = Crop(Deriv(sol.Vr.grid, sol.Phi, 1), [0 1]);
    Et = Crop(Deriv(sol.Vt.grid, sol.Phi, 2), [1 0]) * '1/r';
    gi = sol.P.grid;
    
    Q = Deriv(gi, 'r^2' * Er, 1) * '1/r^2' + ...
        Deriv(gi, 'sin(t)' * Et, 2) * '1/(r * sin(t))';

    gr = Grid(sol.Vr.grid.r(2:end-1), sol.Vr.grid.t(2:end-1));
    gt = Grid(sol.Vr.grid.r(2:end-1), sol.Vt.grid.t);
    
    forceR = Deriv(gr, - sol.P + Deriv(gi, Crop(sol.Vr, [0 1]) * 'r^2', 1) * '1/r^2', 1) ...
        + Deriv(gr, (Deriv(gt, Crop(sol.Vr, [1 0]), 2) - 2*Interp(gt, sol.Vt)) * 'sin(t)', 2) * ('1/(r^2 * sin(t))');
    
    gt = Grid(sol.Vt.grid.r(2:end-1), sol.Vt.grid.t(2:end-1));
    gr = Grid(sol.Vr.grid.r, sol.Vt.grid.t(2:end-1));
    
    forceT = Deriv(gt, sol.P, 2) * '-1/r' ...
        + '1/r^2' * Deriv(gt, 'r^2'*Deriv(gr, Crop(sol.Vt, [0 1]), 1), 1) ...
        + '1/r^2' * Deriv(gt, '1/sin(t)' * Deriv(gi, Crop(sol.Vt, [1 0]) * 'sin(t)', 2) + 2*Interp(gi, sol.Vr), 2);
end

function [eqn] = equations(sol)    
    [forceR, forceT] = momentum(sol);
    eqn = Join(charge(sol), salt(sol), mass(sol), forceR, forceT);
end

function sol = Solution(grid, init)
    sol.grid = grid;    
    c = struct2cell(init);
    sol.x = Variable(c{:});
    sol.numel = numel(sol.x.value);
    offset = 0;
    F = fieldnames(init);
    for k = 1:numel(F)
        name = F{k};
        g = grid.(name); % grid for variable
        m = g.numel; % dimension of variable        
        L = sparse(1:m, offset + (1:m), 1, ...
                     m, sol.numel); % Restrictor
        op = Linear(g, sol.x, L); % Linear operator
        sol.(name) = op;
        offset = offset + g.numel;
    end    
end
