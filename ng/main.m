% Coupled simulation main script.
function main

% Create grids
    nx = 3;
    ny = 3;
    x = linspace(0, 1, nx);
    y = linspace(0, 1, ny);
    [center, interior, xstag, ystag] = grids(x, y);

% Initialize variables
    Phi = zeros(center.sz);
    C = ones(center.sz);
    Vx = zeros(xstag.sz);
    Vy = zeros(ystag.sz);
    P = zeros(interior.sz);
    
    problem = @() struct('iters', 0);
    
    laplace = problem();
    advect = problem();
    stokes = problem();

% Laplace problem (Phi)
    function update_laplace()
        beta = 1;
        gridPhi = center;
        laplace_full_operator = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y, C);

        [lhs1, rhs1] = dirichlet(gridPhi, [-1 0]);
        [lhs2, rhs2] = neumann(gridPhi, [+1 0]);
        [lhs3] = neumann(gridPhi, [0 -1]);
        [lhs4] = neumann(gridPhi, [0 +1]);
        lhs = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridPhi.I);
        laplace.operator = laplace_full_operator * lhs;
        E = -beta*cos(gridPhi.y(2:end-1));
        laplace.rhs = laplace_full_operator * ...
            (rhs1(-log( C(2, 2:end-1) )) + rhs2(E));
        laplace.precond = redblack(interior.sz, laplace.operator);
    end
     
% Diffusion-advection problem (C)
    function update_advection
        alpha = 1;
        gridC = center;
        advection_full_operator = laplacian(gridC.I, gridC.X, gridC.Y) - ...
            advection(gridC.I, gridC.X, gridC.Y, alpha*Vx, alpha*Vy, 'central');

        [lhs1, rhs1] = dirichlet(gridC, [-1 0]);
        [lhs2, rhs2] = dirichlet(gridC, [+1 0]);
        [lhs3] = neumann(gridC, [0 -1]);
        [lhs4] = neumann(gridC, [0 +1]);
        lhs = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridC.I);
        advect.operator = advection_full_operator * lhs;
        advect.rhs = advection_full_operator * ...
            (rhs1(exp( -Phi(2, 2:end-1))) + rhs2(1));
        advect.precond = redblack(interior.sz, advect.operator);
    end    

% Stokes problem
    function update_stokes
        vel = 1;    
        gamma = 0.5;
        L = laplacian(center.I, center.X, center.Y);
        grad = @(dir) gradient(shift(interior.I, -dir), interior.X, interior.Y, dir);
        Gx = grad([1 0]);
        Gy = grad([0 1]);
        q = L * Phi(:);
        Ex = Gx * Phi(center.I);
        Ey = Gy * Phi(center.I);
        
        % P (no boundary condition: 1 DOF, arbirary mean)
        gridP = interior;
        gridVx = xstag;
        gridVy = ystag;
        
        [Gx_p, L_vx, Gx_vx] = stokes1(1, gridVx, gridP); % X dimension
        [Gy_p, L_vy, Gy_vy] = stokes1(2, gridVy, gridP); % Y dimension
        
        [lhs1x, rhs1x] = dirichlet(gridVx, [-1 0]);
        [lhs2x, rhs2x] = dirichlet(gridVx, [+1 0]);
        [lhs3x] = neumann(gridVx, [0 -1]);
        [lhs4x] = neumann(gridVx, [0 +1]);
        lhsVx = lhs1x * lhs2x * lhs3x * lhs4x * restrict(gridVx.I);
        rhsVx =  + rhs1x(0) + rhs2x(-vel*cos(gridVx.y(2:end-1)));
        % (q*Ex)
        [lhs1y, rhs1y] = average_dirichlet(gridVy, [-1 0]);
        [lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
        [lhs3y] = dirichlet(gridVy, [0 -1]);
        [lhs4y] = dirichlet(gridVy, [0 +1]);
        lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
        rhsVy =  + rhs1y(slip(Phi, C, gamma)) + rhs2y(vel*sin(gridVy.y(2:end-1)));
        % (q*Ey)
        div = zeros(gridP.sz); % zero-divergence
        stokes.rhs = [rhsVx(:); rhsVy(:); div(:)];
        
        L_vx1 = L_vx * lhsVx;
        L_vy1 = L_vy * lhsVy;
        Gx_vx1 = Gx_vx * lhsVx;
        Gy_vy1 = Gy_vy * lhsVy;
        stokes.operator = [... % XXX
            L_vx1, sparse(size(L_vx1, 1), size(L_vy1, 2)), -Gx_p; ...
            sparse(size(L_vy1, 1), size(L_vx1, 2)), L_vy1, -Gy_p; ...
            Gx_vx1, Gy_vy1, sparse(prod(gridP.sz), prod(gridP.sz))];
        stokes.precond = stokes_vanka_redblack(gridP.sz, stokes.operator);
    end

% Iterate on problems
    iters = 100;    
    for iter = 1:iters
        update_laplace;         
        Phi = iterate(laplace, Phi);
        
        update_advection;       
        C = iterate(advect, C);

        update_stokes;
        [Vx, Vy, P] = iterate(stokes, Vx, Vy, P);
        P = P - mean(P(:));
    end    
end

function v = slip(varargin)
    v = 0;
end