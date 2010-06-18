% Coupled simulation main script.
function main
    % clc;
% Create grids
    nx = 9;
    ny = 9;
    x = linspace(1, 4, nx);
    y = linspace(0, pi, ny);
    [center, interior, xstag, ystag] = grids(x, y);

% Initialize variables
    Phi = zeros(center.sz);
    C = ones(center.sz);
    Vx = zeros(xstag.sz);
    Vy = zeros(ystag.sz);
    P = zeros(interior.sz);
    
    problem = @(iters) struct('iters', iters);

% Laplace problem (Phi)
    function update_laplace()
        beta = 1;
        gridPhi = center;
        laplace_full_operator = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y, C);

        [lhs1, rhs1] = dirichlet(gridPhi, [-1 0]);
        [lhs2, rhs2] = neumann(gridPhi, [+1 0]);
        [lhs3] = neumann(gridPhi, [0 -1]);
        [lhs4] = neumann(gridPhi, [0 +1]);
        lhsPhi = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridPhi.I);
        laplace.operator = laplace_full_operator * lhsPhi;
        Er = -beta*cos( gridPhi.y(2:end-1) );
        rhsPhi = rhs1(-log( C(2, 2:end-1) )) + rhs2(Er); 
        laplace.rhs = laplace_full_operator * rhsPhi;
        laplace.precond = redblack(interior.sz, laplace.operator);
        laplace.expand = @(Phi) ...
            reshape(lhsPhi * Phi(:) - rhsPhi, gridPhi.sz);
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
        lhsC = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridC.I);
        advect.operator = advection_full_operator * lhsC;
        rhsC = rhs1(exp( -Phi(2, 2:end-1))) + rhs2(1);
        advect.rhs = advection_full_operator * rhsC;
        advect.precond = redblack(interior.sz, advect.operator);
        advect.expand = @(C) ...
            reshape(lhsC * C(:) - rhsC, gridC.sz);
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
        Fx = 0; %(Ax * q) .* Ex;
        Fy = 0; %(Ay * q) .* Ey;
        
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
        rhsVx = rhs1x(0) + rhs2x(-vel*cos(gridVx.y(2:end-1)));

        [lhs1y, rhs1y] = average_dirichlet(gridVy, [-1 0]);
        [lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
        [lhs3y] = dirichlet(gridVy, [0 -1]);
        [lhs4y] = dirichlet(gridVy, [0 +1]);
        lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
        Vt = ddslip(Phi, C, gamma, center.y);
        rhsVy = rhs1y(Vt) + rhs2y(vel*sin(gridVy.y(2:end-1)));
        
        div = zeros(gridP.sz); % zero-divergence
        stokes.rhs = [...
            L_vx * rhsVx(:) + Fx; ...
            L_vy * rhsVy(:) + Fy; ...
            Gx_vx * rhsVx(:) + Gy_vy * rhsVy(:) + div(:)];
        
        L_vx1 = L_vx * lhsVx;
        L_vy1 = L_vy * lhsVy;
        Gx_vx1 = Gx_vx * lhsVx;
        Gy_vy1 = Gy_vy * lhsVy;
        stokes.operator = [... 
            L_vx1, sparse(size(L_vx1, 1), size(L_vy1, 2)), -Gx_p; ...
            sparse(size(L_vy1, 1), size(L_vx1, 2)), L_vy1, -Gy_p; ...
            Gx_vx1, Gy_vy1, sparse(prod(gridP.sz), prod(gridP.sz))];
        stokes.precond = stokes_vanka_redblack(gridP.sz, stokes.operator);
        stokes.expandVx = @(Vx) ...
            reshape(lhsVx * Vx(:) - rhsVx, gridVx.sz);
        stokes.expandVy = @(Vy) ...
            reshape(lhsVy * Vy(:) - rhsVy, gridVy.sz);
    end

% Problems:
    laplace = problem(2);
    advect = problem(2);
    stokes = problem(1);
% Iterate on problems
    iters = 3000;
    for iter = 1:iters
        update_laplace;         
        Phi = iterate(laplace, Phi(2:end-1, 2:end-1));
        Phi = laplace.expand(Phi);
        
        update_advection;       
        C = iterate(advect, C(2:end-1, 2:end-1));
        C = advect.expand(C);

        update_stokes;
        [Vx, Vy, P] = iterate(...
            stokes, Vx(2:end-1, 2:end-1), Vy(2:end-1, 2:end-1), P);
        [Vx] = stokes.expandVx(Vx);
        [Vy] = stokes.expandVy(Vy);
        P = P - mean(P(:));
    end    
    figure(1); clf; show(Phi(:, 2:end-1), 'Phi')
    figure(2); clf; show(C(:, 2:end-1), 'C')
    figure(3); clf; show(Vx(:, 2:end-1), 'V_R')
    figure(4); clf; show(Vy(2:end-1, :), 'V_\Theta')
    figure(5); clf; show(P, 'P')
    save results
end
