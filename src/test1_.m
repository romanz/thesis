function test1_
tic;
    
sz = [20 50];
g = Grid(logspace(0, 1, sz(1)), linspace(0, pi, sz(2)));
% variables with boundary conditions
grid.Phi = g.init('central', 'central');
grid.C   = g.init('central', 'central');
grid.Vr  = g.init('', 'central');
grid.Vt  = g.init('central', '');
grid.P   = g.init('interior', 'interior');

% interior grid (ghost points removed)
igrid = noghost(grid);
sol = Solution(igrid, ...
    zeros(igrid.Phi.sz), ...
    ones(igrid.C.sz), ...
    zeros(igrid.Vr.sz), ...
    zeros(igrid.Vt.sz), ...
    zeros(igrid.P.sz));

[bnd] = boundary(grid, sol);
[eqn] = equations(bnd);
dx = eqn.grad() \ eqn.res();
sol.step(-dx);
end

function sol = boundary(grid, sol)
    Phi = Interp(grid.Phi, sol.Phi);
    C = Interp(grid.C, sol.C);
    Vr = Interp(grid.Vr, sol.Vr);
    Vt = Interp(grid.Vt, sol.Vt);
    P = Interp(grid.P, sol.P);

    dr = diff(grid.Phi.r(end-1:end));
    dPhi = @(r, t) -dr * sol.beta * cos(t);
    
    c0     = Boundary(C,   [-1 0], exp(-sol.Phi('1', ':')));
    phi0   = Boundary(Phi, [-1 0], -log(sol.C('1', ':')));
    
    cInf   = Boundary(C,   [+1 0], 1);
    phiInf = Boundary(Phi, [+1 0], sol.Phi('end', ':') + dPhi);

    VrInf  = Boundary(Vr, [+1, 0], @(r,t) -sol.Vinf*cos(t));
    VtInf  = Boundary(Vt, [+1, 0], @(r,t)  sol.Vinf*sin(t));

    sol.Phi = Symm(Phi + phi0 + phiInf);
    sol.C = Symm(C + c0 + cInf);
    
    g = Grid(1, sol.Vt.grid.t);
    xi = log(Interp(g, sol.C)) - log(sol.gamma);
    
    phi = Interp(Grid(1, sol.Phi.grid.t(2:end-1)), sol.Phi); % Phi at R=1
    Dphi = Deriv(g, phi, 2); % Dphi/Dtheta at R=1
    Vs = Func(xi, '4*log((exp(x/2) + 1)/2)') * Dphi;
    Vt0    = Boundary(Vt, [-1, 0], 2*Vs - sol.Vt('1', ':') );
    
    sol.Vr = Symm(Vr + VrInf);
    sol.Vt = Vt + Vt0 + VtInf;
    sol.P = P;
end

function Boundary(op, dir, val)
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
    
    op = Linear(op.grid, op);
    op.L = sparse([I; J1], [I; J2], 1, op.grid.numel, op.grid.numel);    
end

function flux = charge(sol)
    DPhi_Dr = Deriv(sol.grid.Vr, sol.Phi, 1);
    DPhi_Dt = Deriv(sol.grid.Vt, sol.Phi, 2);
    Cr = Interp(DPhi_Dr.grid, sol.C);
    Ct = Interp(DPhi_Dr.grid, sol.C);
    flux = ...
        Deriv((@(r,t) r^2) * Cr * DPhi_Dr, 1) * @(r,t) r^(-2) + ...
        Deriv((@(r,t) sin(t)) * Ct * DPhi_Dt, 2) * @(r,t) 1/(r^2 * sin(t));
end

function eqn = equations(sol)
    eqn = [charge(sol); salt(sol); force(sol); mass(sol)];
end

function sol = Solution(grid, varargin)
    sol.grid = grid;
    sol.v = Variable(varargin{:});
    sol.numel = numel(sol.x);
    offset = 0;
    F = fieldnames(grid);
    for k = 1:numel(F)
        name = F{k};
        g = grid.(name); % grid for variable
        m = g.numel; % dimension of variable
        G = sparse(1:m, offset + (1:m), 1, m, sol.numel); % Restrictor
        sol.(name) = Linear(g, sol.v, G); % Linear operator
        offset = offset + g.numel;
    end
    function step(dx)
        sol.v.value = sol.v.value + dx;
    end
    sol.step = @step;
end

function grid = noghost(grid)
    for n = fieldnames(grid).'
        g = grid.(n{1}); % remove ghost points from grid
        r = g.r((1+g.ghost(1)):(end-g.ghost(1)));
        t = g.t((1+g.ghost(2)):(end-g.ghost(2)));
        grid.(n{1}) = Grid(r, t);
    end
end
