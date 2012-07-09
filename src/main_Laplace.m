function [sol] = main_Laplace(alpha)
    clc;
    g = grids(1e3, 100);
    
    init.Phi = zeros(g.Phi.size);

    sol = Solution(g, init);
    sol.Du = 0;
    sol.grid = g;
    
    if nargin < 1
        alpha = 0;
    end

    betas = 1;
    tic;
    for k = 1:numel(betas)
        fprintf('==================================================================\n')
        sol.beta = betas(k);
        for i = 1
            sol = solver(sol, [7 0], alpha);
        end
    end
    toc
    r = sol.Phi.grid.r; t = sol.Phi.grid.t; phi = regrid(sol.Phi);
    semilogx(convn(r, [1;1]/2, 'valid'), ...
        diff(phi, 1) ./ repmat(diff(r), [1 size(phi, 2)]), g.r(end), betas(end)*(t(2:end-1)), '.')
end

function sol = solver(sol, iters, alpha)
    [sol, iter] = update(sol); % update equations and iterators
    for k = 1:iters(1)
        
%         for j = 1:iters(2)
%             [r, dx] = iter.bnd();
%             fprintf('\t%e -> %e\n', norm(r), norm(dx))
%         end
        [r, dx] = iter.int();
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
    [sol.bnd, sol.I] = boundary_conditions(sol);
    sol.eqn = system_equations(sol);
    
    iter.int = update_interior(sol);
    iter.bnd = update_boundary(sol);
end

function iter = update_interior(sol)
    var = sol.var;
    bnd = sol.bnd;
    eqn = sol.eqn;
    n = numel(var.value)-1;
    T = sparse(1:n, 1:n, 1, n, n+1); % T = 1;
    function [r, dx] = iter_interior()
        Gb = bnd.grad();
        rb = bnd.res();
        
        Gi = eqn.grad();
        ri = eqn.res();
        
        G = [Gb; Gi];
        r = [rb; ri];
        
        A = T*G*T';
        b = T*r;
        dx = linsolve(A, b);
        dx = T'*dx;
        var.update(-dx);
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
function [g] = grids(Rmax, Nr)
    r = logspace(0, log10(Rmax), Nr);
    dr = diff(r(1:2));
    dt = dr;
    Nt = ceil(2 / dt) + 1;
    Nt = Nt + ~mod(Nt, 2);
    t = linspace(-1, 1, Nt);
    r = r(:);
    t = t(:);
    rg = [r(1)^2 / r(2); r; r(end)^2/r(end-1)];
    tg = [2*t(1) - t(2); t; 2*t(end) - t(end-1)];
    rc = average(rg, [1; 1]/2);
    tc = average(tg, [1; 1]/2);
    
    g.Phi = Grid(rc, tc);
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

function [op, I] = boundary_conditions(sol)
    g = sol.grid.Phi;
    g = Grid(mean(g.r(1:2)), g.t);
    g = g.crop(0, 1);
    phi1 = Boundary(sol.Phi, -1, 2);
    bnd{1} = -Deriv(g, phi1, 1) - @(r,t) 0*sol.beta*t; 

    g = sol.grid.Phi;
    g = Grid(mean(g.r(end-1:end)), g.t);
    g = g.crop(0, 1);
    E = -Deriv(g, Boundary(sol.Phi, 1, 2), 1);

    bnd{2} = E - @(r,t) sol.beta*t; 
    bnd{3} = Symm(sol.Phi, 2);
    
    op = Join(bnd{:});
    I = [interior(sol.Phi, 1)];
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
    g = sol.Phi.grid;
    gr = Grid(sol.grid.r, g.t);
    gt = Grid(g.r, sol.grid.t);
    DPhi_Dr = Crop(Deriv(gr, sol.Phi, 1), [0 1]);
    DPhi_Dt = Crop(Deriv(gt, sol.Phi, 2), [1 0]);
    g = Grid(g.r(2:end-1), g.t(2:end-1));
    fluxR = Deriv(g, 'r^2' * DPhi_Dr, 1);
    fluxT = Deriv(g, '1-t^2' * DPhi_Dt, 2);
    flux = fluxR + fluxT;
end


% Electrokinetical system PDEs
function [eqn] = system_equations(sol)    
    eqn = charge(sol);
end
