function [] = main
    % variables with boundary conditions
    [grid.full, grid.int] = grids(logspace(0, 7, 400), linspace(0, pi, 50));
    
    % interior grid (ghost points removed)
    sol = Solution(grid.int);
    sol.alpha = 0;
    sol.beta = 0.001;
    sol.gamma = 0.1;
    sol.Vinf = 0.001;
    
    iters = 10;
    [full] = boundary(grid, sol);
    [eqn] = equations(full);
    [~, J] = find(sol.P.L); % Pressure indices
    for k = 1:iters
        r = eqn.res();
        norm(r(:), inf)
        G = eqn.grad();
        G = G + 1e-10 * speye(size(G));
        dx = G\r;        
        sol.x.value = sol.x.value - dx;
        sol.x.value(J) = sol.x.value(J) - mean(sol.x.value(J));
    end
end

function [full, interior] = grids(r, t)
    r = r(:);
    t = t(:);
    rg = [2*r(1) - r(2); r; 2*r(end) - r(end-1)];
    tg = [2*t(1) - t(2); t; 2*t(end) - t(end-1)];
    rc = average(rg, [1; 1]/2);
    tc = average(tg, [1; 1]/2);
    
    full.Phi = Grid(rc, tc);
    full.C = Grid(rc, tc);
    full.Vr = Grid(r, tc);
    full.Vt = Grid(rc, t);
    full.P = Grid(rc(2:end-1), tc(2:end-1));
    
    names = fieldnames(full);
    for k = 1:numel(names)
        n = names{k};
        g = full.(n);
        if ~strcmp(n, 'P')
            g = Grid(g.r(2:end-1), g.t(2:end-1));
        end
        interior.(n) = g;
    end
end

function sol = boundary(grid, sol)
    Phi = Interp(grid.Phi, sol.Phi);
    C = Interp(grid.C, sol.C);
    Vr = Interp(grid.Vr, sol.Vr);
    Vt = Interp(grid.Vt, sol.Vt);
    P = Interp(grid.P, sol.P);

    dr = diff(grid.Phi.r(end-1:end));
    dPhi = @(r, t) -dr * sol.beta * cos(t);
    
    c1     = Interp(Grid(sol.C.grid.r(1), sol.C.grid.t), sol.C);
    phi1   = Interp(Grid(sol.Phi.grid.r(1), sol.Phi.grid.t), sol.Phi);
    c0     = Boundary(C,   [-1 0], Func(phi1, 'exp(-phi)'));
    phi0   = Boundary(Phi, [-1 0], Func(c1, '-log(c)'));
    
    cInf   = Boundary(C,   [+1 0], 1);
    phiInf = Interp(Grid(sol.Phi.grid.r(end), sol.Phi.grid.t), sol.Phi);
    phiInf = Boundary(Phi, [+1 0], phiInf + dPhi);

    VrInf  = Boundary(Vr, [+1 0], @(r,t)-sol.Vinf*cos(t));
    VtInf  = Boundary(Vt, [+1 0], @(r,t) sol.Vinf*sin(t));

    sol.Phi = Symm(Phi + phi0 + phiInf);
    sol.C = Symm(C + c0 + cInf);
    
    g = Grid(1, sol.Vt.grid.t);
    xi = -Interp(g, sol.Phi) - log(sol.gamma);
    
    phi   = Interp(Grid(1, sol.Phi.grid.t(2:end-1)), sol.Phi); % Phi at R=1
    Dphi  = Deriv(g, phi, 2); % Dphi/Dtheta at R=1
    Vs    = Func(xi, '4*log((exp(x/2) + 1)/2)') * Dphi;
    Vt1   = Interp(Grid(sol.Vt.grid.r(1), sol.Vt.grid.t), sol.Vt);
    Vt0   = Boundary(Vt, [-1, 0], 2*Vs - Linear(Vs.grid, Vt1, []) );
    
    sol.Vr = Symm(Vr + VrInf);
    sol.Vt = Vt + Vt0 + VtInf;
    sol.P = P;
end

function op = Boundary(op, dir, val)
    g = op.grid; % result grid
    r = F(g.r, dir(1));
    t = F(g.t, dir(2));
    b = Grid(r, t); % boundary grid
    if ~isa(val, 'Operator')
        val = Const(b, val);
    else
        val = Linear(b, val, []);
    end
    op = Interp(g, val);
    function x = F(x, d)
        if d > 0, x = x(end); end
        if d < 0, x = x(1); end
        if d == 0, x = x(2:end-1); end
    end
end

function op = Symm(op)
    I = true(op.grid.sz);
    I(:, [1 end]) = false;
    I = find(I);
    
    J1 = false(op.grid.sz);
    J1(:, [1 end]) = true;
    J1 = find(J1); % ghost points
    
    J2 = false(op.grid.sz);    
    J2(:, [2 end-1]) = true;
    J2 = find(J2); % interior points
    
    L = sparse([I; J1], [I; J2], 1, op.grid.numel, op.grid.numel);
    op = Linear(op.grid, op, L);
end

function flux = charge(sol)
    DPhi_Dr = Crop(Deriv(sol.Vr.grid, sol.Phi, 1), [0 1]);
    DPhi_Dt = Crop(Deriv(sol.Vt.grid, sol.Phi, 2), [1 0]);
    Cr = Interp(DPhi_Dr.grid, sol.C);
    Ct = Interp(DPhi_Dt.grid, sol.C);
    g = sol.Phi.grid;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, (@(r,t) r^2) * Cr * DPhi_Dr, 1) * @(r,t) r^(-2);
    fluxT = Deriv(g, (@(r,t) sin(t)) * Ct * DPhi_Dt, 2) * @(r,t) 1/(r^2 * sin(t));
    flux = fluxR + fluxT;
end

function flux = salt(sol)
    DC_Dr = Crop(Deriv(sol.Vr.grid, sol.C, 1), [0 1]);
    DC_Dt = Crop(Deriv(sol.Vt.grid, sol.C, 2), [1 0]);
    % XXX Add advection!
    g = sol.C.grid;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, (@(r,t) r^2) * DC_Dr, 1) * @(r,t) r^(-2);
    fluxT = Deriv(g, (@(r,t) sin(t)) * DC_Dt, 2) * @(r,t) 1/(r^2 * sin(t));
    flux = fluxR + fluxT;
end

function flux = mass(sol)
    Dm_Dr = Crop(sol.Vr, [0 1]);
    Dm_Dt = Crop(sol.Vt, [1 0]);
    g = sol.P.grid;
    fluxR = Deriv(g, (@(r,t) r^2) * Dm_Dr, 1) * @(r,t) r^(-2);
    fluxT = Deriv(g, (@(r,t) sin(t)) * Dm_Dt, 2) * @(r,t) 1/(r * sin(t));
    flux = fluxR + fluxT;
end

function [forceR, forceT] = momentum(sol)
    Er = Crop(Deriv(sol.Vr.grid, sol.Phi, 1), [0 1]);
    Et = Crop(Deriv(sol.Vt.grid, sol.Phi, 2), [1 0]) * @(r,t) 1/r;
    gi = sol.P.grid;
    
    Q = Deriv(gi, (@(r,t) r^2) * Er, 1) * (@(r,t) r^(-2)) + ...
        Deriv(gi, (@(r,t) sin(t)) * Et, 2) * (@(r,t) 1/(r * sin(t)));

    gr = Grid(sol.Vr.grid.r(2:end-1), sol.Vr.grid.t(2:end-1));
    gt = Grid(sol.Vr.grid.r(2:end-1), sol.Vt.grid.t);
    forceR = Deriv(gr, - sol.P + Deriv(gi, Crop(sol.Vr, [0 1]) * @(r,t)r^2, 1) * @(r,t)1/r^2, 1) ...
        + Deriv(gr, (Deriv(gt, Crop(sol.Vr, [1 0]), 2) - 2*Interp(gt, sol.Vt)) * @(r,t)sin(t), 2) * (@(r,t) 1/(r^2 * sin(t)));
    
    gt = Grid(sol.Vt.grid.r(2:end-1), sol.Vt.grid.t(2:end-1));
    gr = Grid(sol.Vr.grid.r, sol.Vt.grid.t(2:end-1));
    forceT = Deriv(gt, sol.P, 2) * (@(r,t) -1/r) ...
        + (@(r,t) 1/r^2) * Deriv(gt, (@(r,t) r^2)*Deriv(gr, Crop(sol.Vt, [0 1]), 1), 1) ...
        + (@(r,t) 1/r^2) * Deriv(gt, (@(r,t)1/sin(t)) * Deriv(gi, Crop(sol.Vt, [1 0]) * @(r,t)sin(t), 2) + 2*Interp(gi, sol.Vr), 2);
end

function eqn = equations(bnd)    
    [forceR, forceT] = momentum(bnd);
    eqn = Join(charge(bnd), salt(bnd), mass(bnd), forceR, forceT);
end

function sol = Solution(grid, varargin)
    sol.grid = grid;
    sol.x = Variable(varargin{:});
    sol.numel = numel(sol.x.value);
    offset = 0;
    F = fieldnames(grid);
    for k = 1:numel(F)
        name = F{k};
        g = grid.(name); % grid for variable
        m = g.numel; % dimension of variable
        op = Linear(g, sol.x, []); % Linear operator
        J = offset + (1:m); % Colum indices
        op.L = sparse(1:m, J, 1, m, sol.numel); % Restrictor
        sol.(name) = op;
        offset = offset + g.numel;
    end    
end
