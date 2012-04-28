function res = main(init)
    grid = grids(1e7, 300, 60);
    
    V = [0.5:0.1:1.5];
    F = zeros(size(V));
    for i = 1:numel(V)
        init.Phi = zeros(grid.Phi.size);
        init.C = ones(grid.C.size);
        init.Vr = zeros(grid.Vr.size);
        init.Vt = zeros(grid.Vt.size);
        init.P = zeros(grid.P.size);
        
        sol = Solution(grid, init);
        sol.alpha = 0.0;
        sol.beta = 0.001;
        sol.Du = 1;
        sol.zeta = 10;
        sol.Vinf = V(i) * sol.beta * (sol.Du*log(16) + sol.zeta)/(1 + 2*sol.Du);
        force = total_force(sol, grid);

        [sol, iter] = update(sol);
        for k = 1:5
            for j = 1:3
                [r, dx] = iter.bnd();
            fprintf('\t%e -> %e\n', norm(r), norm(dx))
            end
            [r, dx] = iter.int();
            fprintf('>>> %e -> %e\n', norm(r), norm(dx))

            if k == 1
                sol.alpha = 0.5;
                [sol, iter] = update(sol);
            end
        end
        F(i) = force();
    end
    plot(V, F, '.')
    res = stfun(sol, fieldnames(init), @(v) regrid(v));
    save main
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
    [sol.bnd, sol.I] = boundary_conditions(sol);
    sol.eqn = system_equations(sol);
    
    iter.int = update_interior(sol);
    iter.bnd = update_boundary(sol);
end

function iter = update_interior(sol)
    var = sol.var;
    I = sol.I;
    bnd = sol.bnd;
    eqn = sol.eqn;
    Pi = select(I)';
    n = nnz(I)-1;
    T = sparse(1:n, 1:n, 1, n, n+1);
    function [r, dx] = iter_interior()
        G = bnd.grad();
        A = G(:, ~I);
        B = G(:, I);
        H = linsolve(A, B);
        % Hx + y = 0
        
        r = eqn.res();
        G = eqn.grad();
        
        Gi = G(:, I); % interior variables
        Gb = G(:, ~I); % boundary variables
        % Gi x + Gb y = -r
        % H  x +    y =  0
        G = Gi - Gb * H;
        
        A = T*G*T';
        b = T*r;
        dx = linsolve(A, b);
        dx = T'*dx;
        var.update(-Pi*dx);
    end
    iter = @iter_interior;
end

function iter = update_boundary(sol)
    var = sol.var;
    I = sol.I;
    bnd = sol.bnd;
    Pb = select(~I)';
    
    function [r, dx] = iter_boundary()
        r = bnd.res();
        G = bnd.grad();
        dx = linsolve(G*Pb, r);
        var.update(-Pb*dx);
    end

    iter = @iter_boundary;
end

% Solve linear system Ax = B
function x = linsolve(A, B)
    sz = size(A);
    assert(sz(1) == sz(2));
    assert(sz(1) == size(B, 1))
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
        offset = offset + g.numel;
    end    
end

% Create problem grids.
function [g] = grids(Rmax, Nr, Nt)
    r = logspace(0, log10(Rmax), Nr);
    t = linspace(0, pi, Nt + ~mod(Nt, 2));
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

function [op, I] = boundary_conditions(sol)

    g = Grid(sol.Vr.grid.r(1), sol.Vr.grid.t);
    g1 = g.crop(0, 1);
    phi1 = Boundary(sol.Phi, -1, 2);
    c1 = Boundary(sol.C, -1, 2);
    bnd{1} = Deriv(g1, phi1 + log(c1), 1); 

    phi2 = Interp(g, sol.Phi);
    c2 = Interp(g, sol.C);
    bnd{2} = Deriv(g1, c1, 1) - sol.Du * surface_laplacian(phi2 - log(c2), g); 
    
    phi_inf = @(r, t) -sol.beta*r*cos(t);
    bnd{3} = Boundary(sol.Phi, 1, 1) - phi_inf; % Phi, C @ R=inf
    bnd{4} = Boundary(sol.C, 1, 1) - 1; % Phi, C @ R=inf
    
    bnd{5} = Boundary(sol.Vr, 1, 1) - ( @(r,t) -sol.Vinf*cos(t) );
    bnd{6} = Boundary(sol.Vt, 1, 1) - ( @(r,t)  sol.Vinf*sin(t) );
    
    bnd{7} = Symm(sol.Phi, 2);
    bnd{8} = Symm(sol.C, 2);
    bnd{9} = Symm(sol.Vr, 2);
    bnd{10} = Symm(sol.Vt, 1);
        
    g = Grid(1, sol.Vt.grid.t(2:end-1));
    zeta = sol.zeta; % - log(Interp(g, sol.C))
     
    g1 = Grid(1, sol.Phi.grid.t(2:end-1));
    phi1 = Interp(g1, sol.Phi); % Phi at R=1
    logc1 = log(Interp(g1, sol.C)); % log(C) at R=1
    D  = @(f) Deriv(g, f, 2); % Dphi/Dtheta at R=1 on Vt grid.
    Vs = zeta * D(phi1) + (log(16) - zeta) * D(logc1); 
    
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

% Salt flux divergence
function flux = salt(sol)
    % Diffusion    
    DC_Dr = Deriv(sol.Vr.grid, sol.C, 1);
    DC_Dt = Deriv(sol.Vt.grid, sol.C, 2);
    
    % Advection (with upwind discretization)
    Cr = Upwind(sol.C, sol.Vr);
    Ct = Upwind(sol.C, sol.Vt);
    
    % Salt fluxes
    Fr = - DC_Dr         + sol.alpha * Cr * sol.Vr;    
    Ft = - DC_Dt * '1/r' + sol.alpha * Ct * sol.Vt;
    
    % Divergence
    g = sol.C.grid;
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2'    * Crop(Fr, [0 1]), 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Crop(Ft, [1 0]), 2) * '1/(r * sin(t))';
    flux = fluxR + fluxT; % Diffusion
end

% Mass flux divergence
function flux = mass(sol)
    Dm_Dr = Crop(sol.Vr, [0 1]);
    Dm_Dt = Crop(sol.Vt, [1 0]);
    g = sol.P.grid;
    fluxR = Deriv(g, 'r^2' * Dm_Dr, 1) * '1/r^2';
    fluxT = Deriv(g, 'sin(t)' * Dm_Dt, 2) * '1/(r * sin(t))';
    flux = fluxR + fluxT;
end

% Momentum flux divergence
function [forceR, forceT] = momentum(sol)
    Er = Crop(Deriv(sol.Vr.grid, sol.Phi, 1), [0 1]);
    Et = Crop(Deriv(sol.Vt.grid, sol.Phi, 2), [1 0]) * '1/r';
    gi = sol.P.grid;
    
    Q = Deriv(gi, 'r^2' * Er, 1) * '1/r^2' + ...
        Deriv(gi, 'sin(t)' * Et, 2) * '1/(r * sin(t))';

    gr = sol.Vr.grid.crop(1, 1);
    gt = Grid(sol.Vr.grid.r(2:end-1), sol.Vt.grid.t);
    
    forceR = Deriv(gr, - sol.P + Deriv(gi, Crop(sol.Vr, [0 1]) * 'r^2', 1) * '1/r^2', 1) ...
        + Deriv(gr, (Deriv(gt, Crop(sol.Vr, [1 0]), 2) - 2*Interp(gt, sol.Vt)) * 'sin(t)', 2) * ('1/(r^2 * sin(t))');
    forceR = forceR + Selector(gr, Er) * Interp(gr, Q);
    
    gt = sol.Vt.grid.crop(1, 1);
    gr = Grid(sol.Vr.grid.r, sol.Vt.grid.t(2:end-1));
    
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
    P = Interp(g, sol.P);
    
    dVr_dr = Deriv(g, Interp(Grid(sol.Vr.grid.r([1 3]), g.t), sol.Vr), 1);
    dVr_dt = Deriv(g, Selector(Grid(g.r, sol.Vr.grid.t),      sol.Vr), 2);
    
    Vt = Interp(g, sol.Vt);
    dVt_dr = Deriv(g, Selector(Grid(sol.Vt.grid.r(1:2), g.t), sol.Vt), 1);
    
    dPhi_dr = Deriv(g, Interp(Grid(sol.Phi.grid.r(1:2), g.t), sol.Phi), 1);
    dPhi_dt = Deriv(g, Interp(Grid(g.r, sol.Phi.grid.t), sol.Phi), 2);
    dPhi_rdt = dPhi_dt * '1/r';
    
    dFr = -P + 2*dVr_dr + 0.5*(dPhi_dr*dPhi_dr - dPhi_rdt*dPhi_rdt);
    dFt = dVt_dr + (dVr_dt - Vt)*'1/r' + dPhi_dr*dPhi_rdt;
    f = dFr * 'cos(t)' - dFt * 'sin(t)'; % local force
    f = f * 'r*sin(t)'; % axial symmetry
    f = optimize(f);
    
    function F = force() 
        F = simpson(g.t, f.res());
    end
    res = @force;
end
