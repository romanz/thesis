% Electrokinetics Non-Linear Solver.
% Usage:
%     betas = 1:0.25:2;     % beta values for continuation solver.
%     sol0 = [];            % Initial solution (empty for beta = 0).
%     figs = [1 2];         % Figure # for plotting the solutions.
%     [sol, grid, prof] = main(betas, sol0, figs); % Run the solver;
%     profview(0, prof);    % Show the profiling data.
%
% [ 1] 1.000 -> 6.026127e-006
% [ 2] 1.250 -> 2.169751e-002
% [ 3] 1.500 -> 1.648105e-003
% [ 4] 1.750 -> 7.391831e-004
% [ 5] 2.000 -> 5.377388e-004
% [ 6] 2.000 -> 4.732835e-004
% [ 7] 2.000 -> 1.164466e-006
% [ 8] 2.000 -> 6.794101e-010
% [ 9] 2.000 -> 9.394985e-017

function [sol, grid, prof] = main(sol, betas, Vinf, gamma, figs)

    % Create grid and initialize the solver
    grid = grids(logspace(0, 3, 60), linspace(0, pi, 30));
    newton_step = solver(grid);

    if isempty(sol) % Initial solution
        sol.Phi = 0*randn(grid.Phi.sz);
        sol.Vx = 0*randn(grid.Vx.sz);
        sol.Vy = 0*randn(grid.Vy.sz);
        sol.P = 0*randn(grid.P.sz);
        sol.C = 1 + 0*randn(grid.C.sz);
    end
    if isempty(betas) % No iterations are possible if no betas.
        return
    end;
    sol.Vinf = Vinf;
    sol.gamma = gamma;

    k = 1; % iteration index
    fprintf('\n');
    profile('on');
    while true % Newton's method with continuation in beta.
        sol.beta = betas(min(k, numel(betas)));
        [sol, u, f] = newton_step(sol);
        assert( all(sol.C(:) > 0), 'Negative C.' );
        
        e = norm(f, 2) / norm(u, 2); % Residual norm
        fprintf('[%2d] %.3f -> %e\n', k, sol.beta, e);
        k = k + 1;
        if e < 1e-12 && k > numel(betas)
            sol = boundaries(sol, grid);
            break;
        end
    end
    profile('off');
    if numel(figs) >= 1
        figure(figs(1));
        show('121', grid.Phi, sol.Phi, '\Phi');
        show('122', grid.C, sol.C, 'C');
        set(gcf, 'Name', 'Numerical Solution')
    end
    if numel(figs) >= 2
        figure(figs(2));
        X = grid.Phi.X;
        Y = grid.Phi.Y;
        show('121', grid.Phi, sol.beta * (0.25*X.^(-2) - X) .* cos(Y), '\Phi');
        show('122', grid.C, 1 + sol.beta * 0.75*X.^(-2) .* cos(Y), 'C');
        set(gcf, 'Name', 'Analytic Solution')
    end
    prof = profile('info');
end

% 3D mesh plot of the solutions
function show(id, grid, sol, msg)
    subplot(id)
    if ~isempty(grid)
        mesh(grid.X, grid.Y, sol)
    else
        mesh(sol)
    end
    title(msg)
end

% Newton solver iteration for given grid.
function step = solver(grid)

    % Divergence, Gradient and Interpolation
    [D1, G1, I1] = operators(grid.center, 1);
    [D2, G2, I2] = operators(grid.center, 2);
    L = sparse_laplacian(grid.center);
    
    [S, Q, E] = stokes(grid);

    S11 = neumanns(grid.Phi, [1 0], [0 -1], [0 1]);
    S22 = neumanns(grid.C, [0 -1], [0 1]);
    S33 = blkdiag(... % Stokes' part Hessian
        neumanns(grid.Vx, [0 -1], [0 1]), ...
        average_dirichlet(grid.Vy, [-1 0]) * expand(grid.Vy.I), ...
        speye(grid.P.numel, grid.P.numel));
    T = S * S33;
    
    % Application of Newton method for given beta
    step = @solver_step;
    function [sol, u, f] = solver_step(sol)
        % Set boundary conditions
        sol = boundaries(sol, grid); 
        
        % Compute gradients and interpolated values
        I1_C = I1 * sol.C(:);
        I2_C = I2 * sol.C(:);
        G1_Phi = G1 * sol.Phi(:);
        G2_Phi = G2 * sol.Phi(:);
        G1_C = G1 * sol.C(:);
        G2_C = G2 * sol.C(:);
        q = Q * sol.Phi(:); % electric charge (on staggered grid)
        e = E * sol.Phi(:); % electric fields (on staggered grid)

        % f = [div(C grad(Phi)); div(grad(C))] -> 0
        f = [D1 * (I1_C .* G1_Phi) + D2 * (I2_C .* G2_Phi); ...
             D1 * G1_C + D2 * G2_C]; ...
        v = [sol.Vx(:); sol.Vy(:); sol.P(:)];
        g = S * v;
        disp(norm(g, 2) / norm(v, 2))

        H = hessian(sol); % Laplace-Advection Hessian
        
        % Apply Newton step
        dw = -blkdiag(T, H) \ [g;f];
        [dVx, dVy, dP, dPhi, dC] = split(dw, ...
            grid.Vx.sz-2, grid.Vy.sz-2, grid.P.sz, grid.Phi.sz-2, grid.C.sz-2);
        
        sol.Phi(grid.Phi.I) = sol.Phi(grid.Phi.I) + dPhi(:);
        sol.C(grid.C.I) = sol.C(grid.C.I) + dC(:);
        sol.Vx(grid.Vx.I) = sol.Vx(grid.Vx.I) + dVx(:);
        sol.Vy(grid.Vy.I) = sol.Vy(grid.Vy.I) + dVy(:);
        sol.P = sol.P + dP;
        meanP = mean(sol.P(:));
        sol.P = sol.P - meanP;
        
        u = [sol.Phi(grid.Phi.I); sol.C(grid.C.I)];

        % Hessian matrix
        function H = hessian(sol)
            % Phi[0] = -log(C[1])
            S12 = dirichlets(grid.C, [-1 0], -1./sol.C(2, 2:end-1));
            % C[0] = exp(-Phi[1])
            S21 = dirichlets(grid.Phi, [-1 0], -exp(-sol.Phi(2, 2:end-1)));

            L1 = D1 * spdiag(I1_C) * G1 + D2 * spdiag(I2_C) * G2;
            L2 = D1 * spdiag(G1_Phi) * I1 + D2 * spdiag(G2_Phi) * I2;
                        
            H = [L1 * S11 + L2 * S21, L1 * S12 + L2 * S22; ...
                 L * S21, L * S22];      
        end
    end
    % Dirichlet for coupled bounadry conditions
    function S = dirichlets(grid, dir, vals)        
        K = find(~grid.I & shift(grid.I, dir)); % ghost points' equations
        I = grid.I & ~shift(grid.I, -dir); 
        J = find(I(2:end-1, 2:end-1)); % near-ghost interior points
        % vals = derivatives of ghost points w.r.t. interior points
        S = sparse(K, J, vals(:), grid.numel, nnz(grid.I));
    end
    % Neumann boundary conditions
    function S = neumanns(grid, varargin)
        S = expand(grid.I);
        for k = 1:numel(varargin)
            d = varargin{k};
            S = neumann(grid, d) * S;
        end
    end
end

% Set boundary conditions
function [sol] = boundaries(sol, grid)
    % Ghost points for R boundaries
    c0 = exp(-sol.Phi(2, 2:end-1));
    phi0 = -log(sol.C(2, 2:end-1));
    dr = grid.Phi.x(end) - grid.Phi.x(end-1);
    Er = -sol.beta * cos(grid.Phi.y(2:end-1)).';
    phi1 = sol.Phi(end-1, 2:end-1) + Er * dr;
    c1 = 1;        
    V_slip = ddslip(sol, grid).';
    Vx_inf = sol.Vinf * cos(grid.Vx.y(2:end-1)).';
    Vy_inf = -sol.Vinf * sin(grid.Vy.y(2:end-1)).';

    sol.C(1, 2:end-1) = c0;
    sol.C(end, 2:end-1) = c1;
    sol.C(1:end, 1) = sol.C(1:end, 2);
    sol.C(1:end, end) = sol.C(1:end, end-1);
    
    sol.Phi(1, 2:end-1) = phi0;
    sol.Phi(end, 2:end-1) = phi1;
    sol.Phi(1:end, 1) = sol.Phi(1:end, 2);
    sol.Phi(1:end, end) = sol.Phi(1:end, end-1);
    
    sol.Vx(1, 2:end-1) = 0;
    sol.Vx(end, 2:end-1) = Vx_inf;
    sol.Vx(1:end, 1) = sol.Vx(1:end, 2); % Symmetry for V_r
    sol.Vx(1:end, end) = sol.Vx(1:end, end-1); % Symmetry for V_r
    
    sol.Vy(1, 2:end-1) = 2*V_slip - sol.Vy(2, 2:end-1);
    sol.Vy(end, 2:end-1) = Vy_inf;
    sol.Vy(1:end, 1) = 0;
    sol.Vy(1:end, end) = 0;
end

%   Dukhin-Derjaguin slip velocity.
function V = ddslip(sol, grid)
    C = col(mean(sol.C(1:2, 2:end-1))); % average on R=1
    Phi = col(mean(sol.Phi(1:2, 2:end-1))); % average on R=1
    
    xi = log(average(C(:), [1;1]/2) / sol.gamma);
    lnC = log(C); 
    
    theta = grid.center.y(2:end-1); % interior cells' centers
    dtheta = diff(theta(:));
    deriv = @(f) diff(f(:)) ./ dtheta; % Derivation by theta.
    V = xi .* deriv(Phi) + 2 * log(1 - tanh(xi/4).^2) .* deriv(lnC);
    % Tangential velocity component.
end

% Create central and staggered grids for Laplace and Stokes problems.
function [grid] = grids(x, y)
    x = x(:); y = y(:); 
    % extended grid (for ghost points)
    xg = [2*x(1) - x(2); x; 2*x(end) - x(end-1)];
    yg = [2*y(1) - y(2); y; 2*y(end) - y(end-1)];
    % cell-centered coordinates
    xc = average(xg, [1; 1]/2);
    yc = average(yg, [1; 1]/2);
    
    % NDGRID convention
    grid.center = init_grid(xc, yc);        
    grid.Phi = grid.center;
    grid.C = grid.center;
	grid.Vx = init_grid(xg(2:end-1), yc); % V_r grid
    grid.Vy = init_grid(xc, yg(2:end-1)); % V_theta grid
    grid.P = init_grid(xc(2:end-1), yc(2:end-1), false); % Pressure grid
end

% Splits x into seperate variables according to specified sizes.
%   [x1, x2, x3] = split(x, sz1, sz2, sz3);
%   where sz{i} = size(x{i}) and x = [x1(:); x2(:); x3(:)];
function [varargout] = split(x, varargin)
    N = numel(varargin);
    varargout = cell(N);
    offset = 0;
    for k = 1:N
        sz = varargin{k};
        m = prod(sz);
        varargout{k} = reshape(x(offset + (1:m)), sz);
        offset = offset + m;
    end
    assert(offset == numel(x), 'Splitting mismatch.');
end
