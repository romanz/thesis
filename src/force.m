% Electrokinetics Non-Linear Force Solver.
% Usage:
%   betas = 1:0.25:2;     % beta values for continuation solver.
%   sol0 = [];            % Initial solution (empty for beta = 0).
%   velocity = 0.5;       % fluid velocity at infinity.
%   [sol, grid, prof] = main(sol0, betas, velocity); % Run the solver;
%   profview(0, prof);    % Show the profiling data.

function [sol, grid, prof] = force(sol, betas, Vinf, varargin)
    assert(nargin > 0);
    conf = defaults(varargin, ...
        'figures', [], ...
        'version', false, ...
        'profile', false, ...
        'verbose', true, ...
        'radius', logspace(0, 3, 60), ...
        'theta', linspace(0, pi, 45), ...
        'gamma', 10, ...
        'alpha', 0, ...
        'maxres', 1e-9, ...
        'miniters', 1, ...
        'maxiters', 100 ...
    );
    
    if conf.version && exist('git', 'file') == 2 % git wrapper function
        [~, ver] = git('--no-pager', 'log -n1 --format=format:"%h (%ci)"');
        fprintf('\nSolver version %s.\n', ver);
        conf.version = ver;
    end

    % Create grid and initialize the solver
    grid = grids(conf.radius, conf.theta);
    newton_step = newton_solver(grid);
    
    if isempty(sol) % Initialize solution
        sol.Phi = zeros(grid.Phi.sz);
        sol.Vx = zeros(grid.Vx.sz);
        sol.Vy = zeros(grid.Vy.sz);
        sol.P = zeros(grid.P.sz);
        sol.C = 1 + zeros(grid.C.sz);
    else
        assert(isequal(sol.grid.radius, grid.radius));
        assert(isequal(sol.grid.theta, grid.theta));
    end
    sol.alpha = conf.alpha;
    sol.gamma = conf.gamma;
    sol.Vinf0 = -2*log((sol.gamma^(1/4) + sol.gamma^(-1/4)) / ...
        (2*sol.gamma^(1/4)));
    if isempty(betas) % No iterations are possible if no betas.
        return
    end;
    sol.Vinf = Vinf;

    k = 1; % iteration index
    if conf.profile
        profile('on');
    end
    while true % Newton's method with continuation in beta.
        sol.beta = betas(min(k, numel(betas)));
        [sol, res] = newton_step(sol);
        assert( all(sol.C(:) > 0), 'Negative C.' );
        
        res = norm(res, 2) / sqrt(numel(res)); % Residual norm
        if conf.verbose
            fprintf('[%d] B = %e, Residual = %e\n', k, sol.beta, res); 
        end
        
        k = k + 1;
        if (k > conf.miniters) && ...
           (res < conf.maxres || k > conf.maxiters) && k > numel(betas)
            break;
        end
    end
    sol.force = total_force(sol, grid);
    sol.grid = grid;
    sol = streamfunc(sol);
    
    if conf.profile
        profile('off');
    end
    prof = profile('info');
end

function [force] = total_force(sol, grid)
    force = struct();
    radius = grid.radius;
    theta = grid.theta;

    Vx0 = sol.Vx(1:2, :); % radial components
    Vy0 = sol.Vy(1:2, :); % tangential components
    P0 = sol.P(1, :); 
    % grid for computation of (V,P) around R=1.
    grid = grids([2 -1; eye(2)]*radius(1:2), theta);
    stokes_operator = stokes(grid);
    n = numel(grid.P.y); % number of interior cells on R=1.

    % Equations extraction: 
    % * include radial forces and divergence for r<1.
    % * exclude tangential forces and divergence for r>1.
    I = [true(nnz(grid.Vx.I), 1); false(nnz(grid.Vy.I), 1); ...
        repmat([true; false], n, 1)];
    stokes_operator1 = stokes_operator(I, :);

    % Variables extraction:
    % * The unknowns are radial velocities and pressure for r<1.
    J = [col([0 0 0; repmat([1 0 0], n, 1); 0 0 0]'); false(grid.Vy.numel, 1); repmat([1; 0], n, 1)];
    A = stokes_operator1 * expand(J);
    % * All other variables have been computed: 
    %   radial velocities (r>=1), pressures (r>1)
    q = stokes_operator1 * [col([zeros(1, size(Vx0, 2)); Vx0]); ...
                            col([zeros(1, size(Vy0, 2)); Vy0; zeros(1, size(Vy0, 2))]); ...
                            col([0*P0(:) P0(:)]')];
    % Solve A*u + q = 0:
    u = A \ -q;
    % Expand radial velocity and pressure for r<1.
    Vx0 = [u([1, 1:n, n])'; Vx0];
    P0 = [u(n+1:end)'; P0];

    function F = midquad(Fr, Ft)
    % Quadrature (midpoint)
        a = grid.P.y'; % Angle axis (from 0 to pi).
        dF = Fr .* cos(a) - Ft .* sin(a); % Projection on axis of symmetry.
        F = sum( dF' .* (2*pi*sin(grid.P.y)) .* diff(theta));    
    end

    %% Newtonian stress
    % Radial force: -P + 2 dVr/dr
    Fr = -mean(P0) ...
         +2*diff(Vx0([1 3], 2:end-1)) / diff(grid.Vx.x([1 3]));
    % Tangential force: dVt/dr + 1/r {dVr/dt - Vt}
    % NOTE: Vr=0 for r=1.
    Ft =  average(diff(Vy0) / diff(radius(1:2)), [1 1]/2) ...
         -average(Vy0(1:2, :), [1 1; 1 1]/4);

    % Quadrature (midpoint)
    force.newton = midquad(Fr, Ft);    

    %% Maxwell stress
    Phi = sol.Phi(1:2, :);
    dPhi_dr = diff(Phi(:, 2:end-1), 1) / diff(radius(1:2));
    dPhi_dtheta = diff(average(Phi, [1 1; 1 1]/4)) ./ diff(theta.');
    Fr = 1/2 * (dPhi_dr.^2 - dPhi_dtheta.^2);
    Ft = dPhi_dr .* dPhi_dtheta;
    force.maxwell = midquad(Fr, Ft);

    %% Total stress
    force.total = force.newton + force.maxwell;

end

% Create central and staggered grids for Laplace and Stokes problems.
function [grid] = grids(x, y)
    x = x(:); grid.radius = x;
    y = y(:); grid.theta = y;
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
