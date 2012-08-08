function [sol] = main(Rmax, Nr, betas)
    g = grids(Rmax, Nr);
    
    init.Phi = zeros(g.Phi.size);
    init.C = ones(g.C.size);
    init.Vr = zeros(g.Vr.size);
    init.Vt = zeros(g.Vt.size);
    init.P = zeros(g.P.size);
    
    sol = Solution(g, init);
    sol.alpha = 0.5;
    sol.Du = 1;
    sol.zeta = 10;
    force = total_force(sol, g);
    
    betas = betas(:);
    V = zeros(size(betas));
    tic;
    for k = 1:numel(betas)
        fprintf('==================================================================\n')
        sol.beta = betas(k);
        Vinf = sol.beta * (sol.Du*log(16) + sol.zeta)/(1 + 2*sol.Du);
        [iter, v] = secant(Vinf * [0.9, 1.1]);
        iters = 5;
        f = zeros(iters, 1);
        for i = 1:iters
            sol.Vinf = v(i);
            sol = solver(sol, [7 0]);
            f(i) = force();
            v(i+1) = iter(f(i));
            fprintf('------------------------------------------------------------------\n')
            fprintf('B = %e \t V = %e \t F = %e\n', betas(k), v(i), f(i))
            fprintf('------------------------------------------------------------------\n')
        end
        V(k) = v(end);
        fname = datestr(now, 'YYYYmmddhhMMss');
        fprintf('Saving to %s...', fname);
        save(fname)
        fprintf('\n');
    end
    toc
end

function sol = solver(sol, iters)
    [sol, iter] = update(sol); % update equations and iterators
    for k = 1:iters
        [r, dx] = iter();
        fprintf('>>> %e -> %e\n', norm(r, inf), norm(dx, inf))
    end
end

function n = count(op)
    n = 1;
    if isa(op, 'Join')
        for k = 1:numel(op.ops)
            n = n + count(op.ops{k});
        end
        return
    end
    if isa(op, 'Binary')
        n = n + count(op.op1);
        n = n + count(op.op2);
    end
    if isa(op, 'Linear') || isa(op, 'Func')
        n = n + count(op.op);
    end
end

function [sol, iter] = update(sol)
    sol.bnd = boundary_conditions(sol);
    sol.eqn = system_equations(sol);
    
    sol.bnd = optimize(sol.bnd);
    sol.eqn = optimize(sol.eqn);
    
    iter = update_solver(sol);
end

function iter = update_solver(sol)
    var = sol.var;
    bnd = sol.bnd;
    eqn = sol.eqn;
    
    n = var.grid.numel;
    J1 = [sol.Phi.grid.numel, n];
    J2 = [bnd.grid.numel+eqn.ops{1}.grid.numel, n];
    P = speye(n); P(:, J1) = [];
    R = speye(n); R(J2, :) = [];
    
    function [r, dx] = iteration()
        Gb = bnd.grad();
        rb = bnd.res();
        
        Gi = eqn.grad();
        ri = eqn.res();
        
        G = [Gb; Gi];
        r = [rb; ri];
        
        % R G P (-dx) = R r
        A = R*G*P;
        b = R*r;
        dx = linsolve(A, b);
        var.update(-P*dx);
    end
    iter = @iteration;
end

% Solve linear system Ax = B
function x = linsolve(A, B)
    sz = size(A);
    assert(sz(1) == sz(2), 'A must be square.');
    assert(sz(1) == size(B, 1), 'A and B must have same number of rows.')
    x = A \ B;
end

% Create the solution vector
function sol = Solution(grid, init)
    c = struct2cell(init);
    sol.var = Variable(c{:});
    sol.numel = numel(sol.var.value);
    offset = 0;
    F = fieldnames(init);
    for k = 1:numel(F)
        name = F{k};
        g = grid.(name); % grid for variable
        m = g.numel; % dimension of variable        
        L = sparse(1:m, offset + (1:m), 1, ...
                     m, sol.numel); % Restrictor
        op = Linear(g, sol.var, L); % Linear operator
        sol.(name) = op;
        offset = offset + m;
    end    
    sol.grid = grid;
end

% Create problem grids.
function [g] = grids(Rmax, Nr)
    r = logspace(0, log10(Rmax), Nr);
    dr = diff(r(1:2));
    dt = dr;
    Nt = ceil(pi / dt) + 1;
    Nt = Nt + ~mod(Nt, 2);
    t = linspace(0, pi, Nt);
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
    
    dr = diff(r(1:2));
    dt = diff(t(1:2));
    fprintf('Nr:Nt = %d x %d\n', Nr, Nt)
    fprintf('dr:dt = %.4f:%.4f\n', dr, (mean(r(1:2) * dt)))
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

% Ls(f) = Dt(sint * Dt(f))/(r^2 sint)
function res = surface_laplacian(op, g)
    f = Interp(g, op);
    g1 = Grid(g.r, conv(g.t, [1 1]/2, 'valid'));
    dfdt_sint = Deriv(g1, f, 2) * 'sin(t)';
    res = Deriv(g.crop(0, 1), dfdt_sint, 2) * '1/(r^2 * sin(t))';
end

function [op] = boundary_conditions(sol)

    bnd = {};
    
    g = Grid(sol.grid.Vr.r(1), sol.grid.Vr.t);
    g1 = g.crop(0, 1);
    phi1 = Boundary(sol.Phi, -1, 2);
    c1 = Boundary(sol.C, -1, 2);    
    bnd{end+1} = Deriv(g1, c1, 1) + Interp(g1, c1) * Deriv(g1, phi1, 1); 

    phi2 = Interp(g, sol.Phi);
    c2 = Interp(g, sol.C);
    bnd{end+1} = Deriv(g1, c1, 1) - sol.Du * surface_laplacian(phi2 - log(c2), g); 
    
    g = sol.grid.Vr;
    g = Grid(g.r(end), g.t(2:end-1));
    E_inf = @(r,t) sol.beta*cos(t);
    E = -Deriv(g, Boundary(sol.Phi, 1, 2), 1);
    C = Interp(g, Boundary(sol.C, 1, 2));

    bnd{end+1} = E - E_inf; 
    bnd{end+1} = C - 1; % R=inf
    
    bnd{end+1} = Boundary(sol.Vr, 1, 1) - ( @(r,t) -sol.Vinf*cos(t) );
    bnd{end+1} = Boundary(sol.Vt, 1, 1) - ( @(r,t)  sol.Vinf*sin(t) );
    
    bnd{end+1} = Symm(sol.Phi, 2);
    bnd{end+1} = Symm(sol.C, 2);
    bnd{end+1} = Symm(sol.Vr, 2);
    bnd{end+1} = Symm(sol.Vt, 1);
        
    g = Grid(1, sol.grid.Vt.t(2:end-1));
    zeta = sol.zeta - log(Interp(g, sol.C));
     
    g1 = Grid(1, sol.grid.Phi.t(2:end-1));
    phi1 = Interp(g1, sol.Phi); % Phi at R=1
    logc1 = log(Interp(g1, sol.C)); % log(C) at R=1
    D  = @(f) Deriv(g, f, 2); % Dphi/Dtheta at R=1 on Vt grid.
    Vs = zeta * D(phi1) + (log(16) - zeta) * D(logc1); 
    
    bnd{end+1} = Interp(g, sol.Vt) - Vs; % Vt @ R=1
    bnd{end+1} = Boundary(sol.Vr, -1, 1); % Vr @ R=1
%}    
    op = Join(bnd{:});

end

% Apply symmetry boundary conditions (for 2nd dimension).
% n = 1, Dirichlet: f = 0.
% n = 2, Neumann: df/dt = 0.
function eq = Symm(op, n)
    dim = 2;
    
    bnd = Boundary(op, +dim, n); % theta = pi
    if n == 1
        eq1 = bnd;
    else
        g = bnd.grid; g = Grid(g.r, mean(g.t));
        eq1 = Deriv(g, bnd, dim);
    end

    bnd = Boundary(op, -dim, n); % theta = 0
    if n == 1
        eq2 = bnd;
    else
        g = bnd.grid; g = Grid(g.r, mean(g.t));
        eq2 = Deriv(g, bnd, dim);
    end
    
    eq = Join(eq1, eq2);
end

% Charge flux divergence:
%   div(C grad(Phi)) = 0.
function flux = charge(sol)
    DPhi_Dr = Crop(Deriv(sol.grid.Vr, sol.Phi, 1), [0 1]);
    DPhi_Dt = Crop(Deriv(sol.grid.Vt, sol.Phi, 2), [1 0]);
    Cr = Interp(DPhi_Dr.grid, sol.C);
    Ct = Interp(DPhi_Dt.grid, sol.C);
    g = sol.grid.Phi;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2' * Cr * DPhi_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Ct * DPhi_Dt, 2) * '1/(r^2 * sin(t))';
    flux = fluxR + fluxT;
end

% Salt flux divergence
function flux = salt(sol)
    % Diffusion    
    DC_Dr = Deriv(sol.grid.Vr, sol.C, 1);
    DC_Dt = Deriv(sol.grid.Vt, sol.C, 2);
    
    % Advection (with upwind discretization)
    Cr = Upwind(sol.C, sol.Vr);
    Ct = Upwind(sol.C, sol.Vt);
    
    % Salt fluxes
    Fr = - DC_Dr         + sol.alpha * Cr * sol.Vr;    
    Ft = - DC_Dt * '1/r' + sol.alpha * Ct * sol.Vt;
    
    % Divergence
    g = sol.grid.C;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2'    * Crop(Fr, [0 1]), 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Crop(Ft, [1 0]), 2) * '1/(r * sin(t))';
    flux = fluxR + fluxT; % Diffusion
end

% Mass flux divergence
function flux = mass(sol)
    Dm_Dr = Crop(sol.Vr, [0 1]);
    Dm_Dt = Crop(sol.Vt, [1 0]);
    g = sol.grid.P;
    fluxR = Deriv(g, 'r^2' * Dm_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Dm_Dt, 2) * '1/(r * sin(t))';
    flux = fluxR + fluxT;
end

% Momentum flux divergence
function [forceR, forceT] = momentum(sol)
    Er = Crop(Deriv(sol.grid.Vr, sol.Phi, 1), [0 1]);
    Et = Crop(Deriv(sol.grid.Vt, sol.Phi, 2), [1 0]) * '1/r';
    gi = sol.grid.P;
    
    Q = Deriv(gi, 'r^2' * Er, 1) * '1/r^2' + ...
        Deriv(gi, 'sin(t)' * Et, 2) * '1/(r * sin(t))';

    gr = sol.grid.Vr.crop(1, 1);
    gt = Grid(sol.grid.Vr.r(2:end-1), sol.grid.Vt.t);
    
    forceR = Deriv(gr, - sol.P + Deriv(gi, Crop(sol.Vr, [0 1]) * 'r^2', 1) * '1/r^2', 1) ...
        + Deriv(gr, (Deriv(gt, Crop(sol.Vr, [1 0]), 2) - 2*Interp(gt, sol.Vt)) * 'sin(t)', 2) * ('1/(r^2 * sin(t))');
    forceR = forceR + Selector(gr, Er) * Interp(gr, Q);
    
    gt = sol.grid.Vt.crop(1, 1);
    gr = Grid(sol.grid.Vr.r, sol.grid.Vt.t(2:end-1));
    
    forceT = Deriv(gt, sol.P, 2) * '-1/r' ...
        + '1/r^2' * Deriv(gt, 'r^2'*Deriv(gr, Crop(sol.Vt, [0 1]), 1), 1) ...
        + '1/r^2' * Deriv(gt, '1/sin(t)' * Deriv(gi, Crop(sol.Vt, [1 0]) * 'sin(t)', 2) + 2*Interp(gi, sol.Vr), 2);
    forceT = forceT + Selector(gt, Et) * Interp(gt, Q);
end

% Electrokinetical system PDEs
function [eqn] = system_equations(sol)    
    [forceR, forceT] = momentum(sol);
    eqn = Join(charge(sol), salt(sol), forceR, forceT, mass(sol));
end

function [res] = total_force(sol, grid)
    % Evaluate total force on the particle by numeric quadrature on r=R.
    g = Grid(grid.r(2), grid.t);
    
    P = Interp(Grid(g.r, sol.grid.Phi.t), sol.P); % interpolate on r=R.
    n = P.grid.numel;
    P = Linear(P.grid, P, sparse(1:n, [2 2:n-1 n-1], 1, n, n));
    P = Interp(g, P);
    
    dVr_dr = Deriv(g, Interp(Grid(sol.grid.Vr.r([1 3]), g.t), sol.Vr), 1);
    dVr_dt = Deriv(g, Selector(Grid(g.r, sol.grid.Vr.t),      sol.Vr), 2);
    
    Vt = Interp(g, sol.Vt);
    dVt_dr = Deriv(g, Selector(Grid(sol.grid.Vt.r(2:3), g.t), sol.Vt), 1);
    
    dPhi_dr = Deriv(g, Interp(Grid(sol.grid.Phi.r(2:3), g.t), sol.Phi), 1);
    dPhi_dt = Deriv(g, Interp(Grid(g.r, sol.grid.Phi.t), sol.Phi), 2);
    dPhi_rdt = dPhi_dt * '1/r';
    
    dFr = -P + 2*dVr_dr + 0.5*(dPhi_dr*dPhi_dr - dPhi_rdt*dPhi_rdt);
    dFt = dVt_dr + (dVr_dt - Vt)*'1/r' + dPhi_dr*dPhi_rdt;
    f = dFr * 'cos(t)' - dFt * 'sin(t)'; % local force
    f = f * 'r*sin(t)'; % axial symmetry
    
    function F = force() 
        F = simpson(g.t, f.res());
    end
    res = @force;
end
