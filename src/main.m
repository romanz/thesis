% Electrokinetics Non-Linear Solver.
% Usage:
%     betas = 1:0.25:2;     % beta values for continuation solver.
%     sol0 = [];            % Initial solution (empty for beta = 0).
%     figs = [1 2];         % Figure # for plotting the solutions.
%     [sol, grid, prof] = main(betas, sol0, figs); % Run the solver;
%     profview(0, prof);    % Show the profiling data.
function [sol, grid, prof] = main(betas, sol, figs)

    % Create grid and initialize the solver
    grid = grids(logspace(0, 3, 100), linspace(0, pi, 60));
    newton_step = solver(grid);

    if isempty(sol) % Initial solution
        sol.Phi = 0*randn(grid.Phi.sz);
        sol.C = 1 + 0*randn(grid.C.sz);
    end
    if isempty(betas) % No iterations are possible if no betas.
        return
    end;

    k = 1; % iteration index
    fprintf('\n');
    profile('on');
    while true % Newton's method with continuation in beta.
        b = betas(min(k, numel(betas)));
        [sol, u, f] = newton_step(sol, b);
        assert( all(sol.C(:) > 0), 'Negative C.' );
        
        e = norm(f, 2) / norm(u, 2); % Residual norm
        fprintf('[%2d] %.3f -> %e\n', k, b, e);
        k = k + 1;
        if e < 1e-12 && k >= numel(betas)
            sol = boundaries(sol, grid, b);
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
        show('121', grid.Phi, b * (0.25*X.^(-2) - X) .* cos(Y), '\Phi');
        show('122', grid.C, 1 + b * 0.75*X.^(-2) .* cos(Y), 'C');
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
    L = D1 * G1 + D2 * G2;
    
    % Application of Newton method for given beta
    step = @solver_step;
    function [sol, u, f] = solver_step(sol, beta)
        % Set boundary conditions
        sol = boundaries(sol, grid, beta); 
        
        % Compute gradients and interpolated values
        I1_C = I1 * sol.C(:);
        I2_C = I2 * sol.C(:);
        G1_Phi = G1 * sol.Phi(:);
        G2_Phi = G2 * sol.Phi(:);
        G1_C = G1 * sol.C(:);
        G2_C = G2 * sol.C(:);
                
        % f = [div(C grad(Phi)); div(grad(C))] -> 0
        f = [D1 * (I1_C .* G1_Phi) + D2 * (I2_C .* G2_Phi); ...
             D1 * G1_C + D2 * G2_C];
        H = hessian(sol);
        du = -(H \ f);
        
        % Apply Newton step on the interior
        [dPhi, dC] = split(du, grid.Phi.sz-2, grid.C.sz-2);
        sol.Phi(grid.Phi.I) = sol.Phi(grid.Phi.I) + dPhi(:);
        sol.C(grid.C.I) = sol.C(grid.C.I) + dC(:);
        u = [sol.Phi(grid.Phi.I); sol.C(grid.C.I)];
        
        % Hessian matrix
        function H = hessian(sol)
            S11 = neumanns(grid.Phi, [1 0], [0 -1], [0 1]);
            S22 = neumanns(grid.C, [0 -1], [0 1]);
            % Phi[0] = -log(C[1])
            S12 = dirichlets(grid.C, [-1 0], -1./sol.C(2, 2:end-1));
            % C[0] = exp(-Phi[1])
            S21 = dirichlets(grid.Phi, [-1 0], -exp(-sol.Phi(2, 2:end-1)));

            L1 = D1 * spdiag(I1_C) * G1 + D2 * spdiag(I2_C) * G2;
            L2 = D1 * spdiag(G1_Phi) * I1 + D2 * spdiag(G2_Phi) * I2;
            H = [L1 * S11 + L2 * S21, ...
                 L1 * S12 + L2 * S22; ...
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
function [sol] = boundaries(sol, grid, beta)
    % Ghost points for R boundaries
    c0 = exp(-sol.Phi(2, 2:end-1));
    phi0 = -log(sol.C(2, 2:end-1));
    dr = grid.Phi.x(end) - grid.Phi.x(end-1);
    Er = -beta * cos(grid.Phi.y(2:end-1)).';
    phi1 = sol.Phi(end-1, 2:end-1) + Er * dr;
    c1 = 1;        
    % Set R boundaries
    sol.C(1, 2:end-1) = c0;
    sol.C(end, 2:end-1) = c1;
    sol.Phi(1, 2:end-1) = phi0;
    sol.Phi(end, 2:end-1) = phi1;
    % Set Theta boundaries
    sol.C(1:end, 1) = sol.C(1:end, 2);
    sol.C(1:end, end) = sol.C(1:end, end-1);
    sol.Phi(1:end, 1) = sol.Phi(1:end, 2);
    sol.Phi(1:end, end) = sol.Phi(1:end, end-1);
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
    grid.P = init_grid(xc(2:end-1), yc(2:end-1)); % Pressure grid
    grid.P.I = true(grid.P.sz); % only interior points are considered
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
    assert(offset == numel(x));
end
