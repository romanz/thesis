% load
% v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
% b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
% loglog(b, v)
function solver1
    tic;
    clc;
    sol = struct( ...
        'radius', logspace(0, 7, 500), ...
        'theta', linspace(0, pi, 50), ...
        'alpha', 1e-3 ...
    );

    sol.grid = grids(sol.radius, sol.theta);
    sol.newton_step = newton_solver1(sol.grid);
    
    sol.Phi = zeros(sol.grid.Phi.sz);
    sol.Vx = cos(sol.grid.Vx.Y);
    sol.Vy = -sin(sol.grid.Vy.Y);
    sol.P = zeros(sol.grid.P.sz);
    sol.C = ones(sol.grid.C.sz);
    
    for k = 1:5
        sol = sol.newton_step(sol);
%         assert( ~any(sol.C(:) < 0), 'C < 0 detected!' );
        
        res = norm(sol.res, inf); % Residual norm
        disp(res)
    end
    g = sol.grid.C;
    c = max(sol.C, 0);
    clf;
    loglog(g.x, c(:, [1 end]), g.x, [exp(-1*g.x*sol.alpha)./g.x]); grid on
%     mesh(x(I, :), y(I, :), (sol.C(I, :)))
%     axis equal
    ylim([1e-10 1])
    title('Diffusion-Advection: -\nabla^2 C + \alpha V \cdot \nabla C = 0')
    xlabel Radius
    ylabel C
    legend('\theta = 0', '\theta = \pi', 'asymptotic')

    save solver1
end

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
