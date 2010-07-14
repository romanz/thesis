% Coupled simulation main script.
function main(filename)
    % clc;
% Create grids
    nx = 17;
    ny = 17;
    x = linspace(1, 4, nx);
    y = linspace(0, pi, ny);
    [center, interior, xstag, ystag] = grids(x, y);

% Initialize variables
    if nargin == 0
        Phi = zeros(center.sz);
        C = ones(center.sz);
        Vx = zeros(xstag.sz);
        Vy = zeros(ystag.sz);
        P = zeros(interior.sz);
    else
        load(filename, 'Phi', 'C', 'Vx', 'Vy', 'P')
    end
    
    problem = @(iters) struct('iters', iters);
% Problems:
    laplace = problem(1);
    advect = problem(1);
    stokes = problem(1);
    iters = 10000;

% Laplace problem (Phi)
    beta = 1;
    gridPhi = center;

    [lhs1, rhs1] = dirichlet(gridPhi, [-1 0]);
    [lhs2, rhs2] = neumann(gridPhi, [+1 0]);
    [lhs3] = neumann(gridPhi, [0 -1]);
    [lhs4] = neumann(gridPhi, [0 +1]);
    lhsPhi = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridPhi.I);
    Er = -beta*cos( gridPhi.y(2:end-1) );
    function update_laplace()
        laplace_full_operator = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y, C); % Spherical
        laplace.operator = laplace_full_operator * lhsPhi;
        rhsPhi = rhs1(-log( C(2, 2:end-1) )) + rhs2(Er); 
        laplace.rhs = laplace_full_operator * rhsPhi;
        laplace.precond = redblack(interior.sz, laplace.operator);
        laplace.expand = @(Phi) ...
            reshape(lhsPhi * Phi(:) - rhsPhi, gridPhi.sz);
    end
     
% Diffusion-advection problem (C)
    alpha = 1;
    gridC = center;
    [lhs1, rhs1] = dirichlet(gridC, [-1 0]);
    [lhs2, rhs2] = dirichlet(gridC, [+1 0]);
    [lhs3] = neumann(gridC, [0 -1]);
    [lhs4] = neumann(gridC, [0 +1]);
    lhsC = (lhs1 * lhs2 * lhs3 * lhs4) * restrict(gridC.I);
    L = laplacian(center.I, center.X, center.Y); % Spherical
    function update_advection
        adv1 = advection(gridC.I, gridC.X, gridC.Y, Vx, Vy, 'upwind'); % Spherical
        advection_full_operator = L - alpha*adv1;
            

        advect.operator = advection_full_operator * lhsC;
        rhsC = rhs1(exp( -Phi(2, 2:end-1))) + rhs2(1);
        advect.rhs = advection_full_operator * rhsC;
        advect.precond = redblack(interior.sz, advect.operator);
        advect.expand = @(C) ...
            reshape(lhsC * C(:) - rhsC, gridC.sz);
    end    

% Stokes problem
    vel = 1;
    gamma = 0.5;
    grad = @(dir) gradient(shift(interior.I, -dir), interior.X, interior.Y, dir); % Spherical
    Gx = grad([1 0]); % at xstag interior
    Gy = grad([0 1]); % at ystag interior
    Ax = interpolator(center.X(2:end-1, 2:end-1), ...
        xstag.X(2:end-1, 2:end-1));
    Ay = interpolator(center.Y(2:end-1, 2:end-1), ...
        ystag.Y(2:end-1, 2:end-1));

    gridP = interior;
    gridVx = xstag;
    gridVy = ystag;

    % Spherical
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

    L_vx1 = L_vx * lhsVx;
    L_vy1 = L_vy * lhsVy;
    Gx_vx1 = Gx_vx * lhsVx;
    Gy_vy1 = Gy_vy * lhsVy;
    
    stokes.operator = [... 
        L_vx1, sparse(size(L_vx1, 1), size(L_vy1, 2)), -Gx_p; ...
        sparse(size(L_vy1, 1), size(L_vx1, 2)), L_vy1, -Gy_p; ...
        Gx_vx1, Gy_vy1, sparse(prod(gridP.sz), prod(gridP.sz))];
    stokes.precond = stokes_vanka_redblack(gridP.sz, stokes.operator);
    
    function update_stokes
        q = L * Phi(:);
        Ex = Gx * Phi(center.I);
        Ey = Gy * Phi(center.I);
        Fx = (Ax * q) .* Ex(:);
        Fy = (Ay * q) .* Ey(:);
        
        % P (no boundary condition: 1 DOF, arbirary mean)
        Vt = ddslip(Phi, C, gamma, center.y);
        rhsVy = rhs1y(Vt) + rhs2y(vel*sin(gridVy.y(2:end-1)));
        
        stokes.rhs = [...
            L_vx * rhsVx(:) - Fx; ...
            L_vy * rhsVy(:) - Fy; ...
            Gx_vx * rhsVx(:) + Gy_vy * rhsVy(:)];
        
        stokes.expandVx = @(Vx) ...
            reshape(lhsVx * Vx(:) - rhsVx, gridVx.sz);
        stokes.expandVy = @(Vy) ...
            reshape(lhsVy * Vy(:) - rhsVy, gridVy.sz);
    end

    dPhi = zeros(iters, 1);
    dC = zeros(iters, 1);
    dVx = zeros(iters, 1);
    dVy = zeros(iters, 1);
    dP = zeros(iters, 1);
    err = @(x) norm(x(:), inf);
    tic;

    % Iterate on problems
    for iter = 1:iters
        update_laplace;         
        Phi1 = iterate(laplace, Phi(2:end-1, 2:end-1));
        Phi1 = laplace.expand(Phi1);
        dPhi(iter) = err(Phi - Phi1);
        Phi = Phi1;
        
        update_advection;       
        C1 = iterate(advect, C(2:end-1, 2:end-1));
        C1 = advect.expand(C1);
        dC(iter) = err(C - C1);
        C = C1;

        update_stokes;
        [Vx1, Vy1, P1] = iterate(...
            stokes, Vx(2:end-1, 2:end-1), Vy(2:end-1, 2:end-1), P);
        [Vx1] = stokes.expandVx(Vx1);
        [Vy1] = stokes.expandVy(Vy1);
        P1 = P1 - mean(P1(:));
        dVx(iter) = err(Vx - Vx1);
        dVy(iter) = err(Vy - Vy1);
        dP(iter) = err(P - P1);
        Vx = Vx1;
        Vy = Vy1;
        P = P1;
    end    
    toc;
    save results
%     figure(1); clf; show(Phi(:, 2:end-1), '\Phi')
%     figure(2); clf; show(C(:, 2:end-1), 'C')
%     figure(3); clf; show(Vx(:, 2:end-1), 'V_R')
%     figure(4); clf; show(Vy(:, :), 'V_\theta')
%     figure(5); clf; show(P, 'P')
%     figure(6); clf; semilogy([dPhi, dC, dVx, dVy, dP]);
%     legend('\Phi', 'C', 'V_R', 'V_\theta', 'P')
end
