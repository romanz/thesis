function F = main1(filename, do_init, iters, beta, gamma, Vinf)
    tic;
    radius = logspace(0, 2, 60).';
    theta = linspace(0, pi, 15).';
    [center, interior, xstag, ystag] = ...
        grids(radius, theta);
    gridPhi = center;
    gridVx = xstag;
    gridVy = ystag;
    gridP = interior;
    gridC = center;

    alpha = 1;    
   
    relax_maxwell = init_maxwell();
    relax_stokes = init_stokes();
    relax_advection = init_advection();    

    if do_init
        solPhi = 0*randn(gridPhi.sz);
        solVx = 0*randn(gridVx.sz);
        solVy = 0*randn(gridVy.sz);
        solP = 0*randn(gridP.sz);
        solC = 1+0*randn(gridC.sz);
    else
        S = load(filename);
        solPhi = S.solPhi;
        solVx = S.solVx;
        solVy = S.solVy;
        solP = S.solP;
        solC = S.solC;
    end
    fprintf('Initialization done after %.3fs.\n', toc);
    
    function state = iteration(state)
        [solPhi, solVx, solVy, solP, solC] = split(state, ...
            gridPhi.sz, gridVx.sz, gridVy.sz, gridP.sz, gridC.sz);
        relax_maxwell(1);
        relax_stokes(1);
        relax_advection(1);        
        state = [solPhi(:); solVx(:); solVy(:); solP(:); solC(:)];
    end

    state = [solPhi(:); solVx(:); solVy(:); solP(:); solC(:)];
    [state] = extrapolate(state, @iteration, 0, 20, 'N/A');
%     beeper(440, 0.1);
    batch = 1;
    [state, err] = extrapolate(state, @iteration, ...
        batch-1, ceil(iters/batch), 'N/A');
    [solPhi, solVx, solVy, solP, solC] = split(state, ...
        gridPhi.sz, gridVx.sz, gridVy.sz, gridP.sz, gridC.sz);
    fprintf('Iterations done after %.3fs.\n', toc);
%     beeper(440, 0.1);
    
    F = total_stress(solVx, solVy, solP, radius, theta);
    fprintf('{%.5e}, {%.5e}\n', F, 6*pi*Vinf);
    save(filename);
    semilogy(err);
% plot([average(-solPhi(1:2, :)', [1 1]/2) ...
%       average(solC(1:2, :)' - 1, [1 1]/2)]);
% plot()

    function func = init_maxwell()
        [P, Q] = maxwell_boundary_cond(gridPhi, beta);
        func = @relax_maxwell;
        function e = relax_maxwell(iters)
            maxwell_operator = maxwell_op(gridPhi, solC); % 
            A = maxwell_operator * P;
            M = redblack(gridPhi.sz-2, A);
            q = Q * maxwell_boundary_vec(solC);
            b = maxwell_operator * q;        
            u = solPhi(gridPhi.I);
            [u, e] = relax(M, A, -b, u, iters);
            solPhi = reshape(P*u + q, gridPhi.sz);
        end
    end

    function func = init_stokes()
        L = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y);
        Lx = interpolator(gridPhi.X(2:end-1, 2:end-1), gridVx.X(2:end-1, 2:end-1)) * L;
        Ly = interpolator(gridPhi.Y(2:end-1, 2:end-1), gridVy.Y(2:end-1, 2:end-1)) * L;
        grad = @(G, dir) gradient(G.I & shift(G.I, -dir), G.X, G.Y, dir);
        Gx = grad(gridPhi, [1 0]);
        Gy = grad(gridPhi, [0 1]);
        divergence = zeros(gridP.numel, 1);
        function rhs = stokes_rhs(solPhi)
            qx = Lx * solPhi(:);
            qy = Ly * solPhi(:);
            Ex = Gx * solPhi(:);
            Ey = Gy * solPhi(:);
            rhs = [qx .* Ex; qy .* Ey; divergence];
        end

        stokes_operator = spheric_stokes(gridVx, gridVy, gridP);
        [P, Q] = stokes_boundary_cond(gridVx, gridVy, gridP, Vinf);
        A = stokes_operator * P;        
        M = stokes_vanka_redblack(gridP.sz, A);
        func = @relax_stokes;
        function e = relax_stokes(iters)
            q = Q * stokes_boundary_vec(solPhi, solC, interior.y, gamma);
            b = stokes_operator * q;        
            u = [solVx(gridVx.I); solVy(gridVy.I); solP(:)];
            [u, e] = relax(M, A, stokes_rhs(solPhi)-b, u, iters);
            u = P*u + q;
            [solVx, solVy, solP] = split(u, gridVx.sz, gridVy.sz, gridP.sz);
            solP = solP - mean(solP(:));
        end
    end


    function func = init_advection()
        [P, Q] = advection_boundary_cond(gridC, 1);
        L = laplacian(gridC.I, gridC.X, gridC.Y);
        func = @relax_advection;
        function e = relax_advection(iters)
            VG = advection(gridC.I, gridC.X, gridC.Y, solVx, solVy);
            advect_operator = L - alpha*VG;
            A = advect_operator * P;
            q = Q * advection_boundary_vec(solPhi);
            b = advect_operator * q;
            M = redblack(gridC.sz-2, A);
            u = solC(gridC.I);
            [u, e] = relax(M, A, -b, u, iters);
            solC = reshape(P*u + q, gridC.sz);
        end
    end
end

function [u, e] = relax(M, A, f, u, iters)
    u0 = u;
    for iter = 1:iters
        for i = 1:numel(M)
            res = f - A*u;
            u = u + M{i}*res;
        end
    end
    e = u - u0;
end

function operator = maxwell_op(gridPhi, solC)
    operator = laplacian(gridPhi.I, gridPhi.X, gridPhi.Y, solC);
end

function [P, Q] = maxwell_boundary_cond(gridPhi, beta)
    [Q1] = dirichlet(gridPhi, [-1 0]); %C
    [P, Q2] = neumann(gridPhi, [+1 0]); %E
    R = symmetry(gridPhi);
    P = R * P * expand(gridPhi.I);
    Einf = -beta*cos(gridPhi.y(2:end-1));
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
