function main1
    clc
    maxwell = struct(); 
    stokes = struct();
    advection = struct();
    
    [center, interior, xstag, ystag] = ...
        grids(logspace(0, .6, 9), linspace(0, pi, 4));
    gridPhi = center;
    gridVx = xstag;
    gridVy = ystag;
    gridP = interior;
    gridC = center;
    
    solPhi = 0.001*randn(gridPhi.sz);
    solVx = 0.001*randn(gridVx.sz);
    solVy = 0.001*randn(gridVy.sz);
    solP = 0.001*randn(gridP.sz);
    solC = 1+0.001*randn(gridC.sz);
    
    alpha = 1;
    beta = 0e-3;
    gamma = 2;
    Vinf = 0e-3;
    iters = 4000;
    for iter = 1:iters
        do_maxwell();
        do_stokes();
        do_advection();
    end
    solPhi, solVx, solVy, solP, solC
    save results

    function do_maxwell()
        maxwell.operator = maxwell_op(gridPhi, solC);
        [P, Q] = maxwell_boundary_cond(gridPhi, beta);
        A = maxwell.operator * P;
        q = Q * maxwell_boundary_vec(solC);
        b = maxwell.operator * q;
        
        M = jacobi(gridPhi.sz-2, A);
        u = solPhi(gridPhi.I);
        u = relax(M, A, -b, u);
        solPhi = reshape(P*u + q, gridPhi.sz);
    end

    function do_stokes()        
        stokes.operator = spheric_stokes(gridVx, gridVy, gridP);
        [P, Q] = stokes_boundary_cond(gridVx, gridVy, gridP, Vinf);
        A = stokes.operator * P;
        q = Q * stokes_boundary_vec(solPhi, solC, interior.y, gamma);
        b = stokes.operator * q;
        
        M = stokes_vanka_redblack(gridP.sz, A);
        u = [solVx(gridVx.I); solVy(gridVy.I); solP(:)];
        u = relax(M, A, -b, u); % TODO: maxwell_forces(Phi)
        u = P*u + q;
        [solVx, solVy, solP] = split(u, gridVx.sz, gridVy.sz, gridP.sz);
        solP = solP - mean(solP(:));
    end

    function do_advection()
        advection.operator = advection_op(gridC, alpha*solVx, alpha*solVy);
        [P, Q] = advection_boundary_cond(gridC, 1);
        A = advection.operator * P;
        q = Q * advection_boundary_vec(solPhi);
        b = advection.operator * q;
        
        M = jacobi(gridC.sz-2, A);
        u = solC(gridC.I);
        u = relax(M, A, -b, u);
        solC = reshape(P*u + q, gridC.sz);
    end
end

function u = relax(M, A, f, u)
    for i = 1:numel(M)
        u = u + M{i}*(f - A*u);
    end
end

function operator = maxwell_op(gridPhi, C)
    operator = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y, C);
end

function operator = advection_op(gridC, Vx, Vy)
    VG = advection(gridC.I, gridC.X, gridC.Y, Vx, Vy, 'central'); %/upwind
    operator = laplacian(gridC.I, gridC.X, gridC.Y) - VG;
end

function [P, Q] = maxwell_boundary_cond(gridPhi, beta)
    [Q1] = dirichlet(gridPhi, [-1 0]); %C
    [P, Q2] = neumann(gridPhi, [+1 0]); %E
    R = symmetry(gridPhi);
    P = R * P * expand(gridPhi.I);
    Einf = beta*cos(gridPhi.y(2:end-1));
    Q = R * [Q1, Q2*Einf];
end
function vec = maxwell_boundary_vec(C)
    vec = [-log(C(2, 2:end-1)).'; 1];
end

function [P, Q] = advection_boundary_cond(gridC, Cinf)
    [Q1] = dirichlet(gridC, [-1 0]); %Phi
    [Q2] = dirichlet(gridC, [+1 0]); %Cinf
    R = symmetry(gridC);
    P = R * expand(gridC.I);
    Cinf = repmat(Cinf, size(Q2, 2), 1);
    Q = R * [Q1, Q2*Cinf];
end
function vec = advection_boundary_vec(Phi)
    C1 = exp(-Phi(2, 2:end-1)).';
    vec = [C1; 1];
end

function [P, Q] = stokes_boundary_cond(gridVx, gridVy, gridP, Vinf)
    R = symmetry(gridVx);
    Px = R * expand(gridVx.I);    
    qx = R * dirichlet(gridVx, [+1 0]) * (-Vinf) * cos(gridVx.y(2:end-1));

    [Py, Qs] = average_dirichlet(gridVy, [-1 0]); % Vslip
    Py = Py * expand(gridVy.I);
    qy = Vinf * dirichlet(gridVy, [+1 0]) * sin(gridVy.y(2:end-1));
    
    P = blkdiag(Px, Py, speye(gridP.numel));
    m = size(Qs, 2);
    Q = [[sparse(gridVx.numel, m); Qs; sparse(gridP.numel, m)], ...
         [qx; qy; sparse(gridP.numel, 1)]];
end
function vec = stokes_boundary_vec(Phi, C, theta, gamma)
    vec = [ddslip(Phi, C, gamma, theta); 1];
end
function R = symmetry(G)
    R = neumann(G, [0 -1], 1) * neumann(G, [0 +1], 1);
end