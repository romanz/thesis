function test

    profile on
    [central] = grids(logspace(0, 2, 90), linspace(0, pi, 50));
    grid.Phi = central;
    grid.C = central;

    newton_step = solver(grid);

    sol.Phi = zeros(grid.Phi.sz);
    sol.C = ones(grid.C.sz);
    betas = [...
        linspace(0.5, 1, 6), ...
        linspace(1.1, 1.3, 11), ...
        linspace(1.3, 1.5, 21)];
%     betas = 0.5;
    k = 1;
    while true
        b = betas(min(k, numel(betas)));
        [sol, u, f] = newton_step(sol, b);
        assert( all(sol.C(:) > 0) );
        e = norm(f, 2) / norm(u, 2);
        fprintf('[%2d] %.3f -> %e\n', k, b, e);
        k = k + 1;
        if e < 1e-10 && k >= numel(betas)
            break;
        end
    end
    profile off
    show(1, grid.Phi, sol.Phi, '\Phi');
    show(2, grid.C, sol.C, 'C');
    save results
    profile viewer
end

function show(id, grid, sol, msg)
    figure(id)
    if ~isempty(grid)
        mesh(grid.X, grid.Y, sol)
    else
        mesh(sol)
    end
    title(msg)
end

function step = solver(grid)
    % Operator definition
    L.Phi = grid_laplacian(grid.Phi);
    L.C = grid_laplacian(grid.C);
    function [sol] = boundaries(sol, beta)
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
    % Newton step
    function [sol, u, f] = step_func(sol, beta)
        sol = boundaries(sol, beta);
        f = [L.Phi * sol.Phi(:); L.C * sol.C(:)];
        H = hessian(sol);
        du = H \ f;
        u = [sol.Phi(grid.Phi.I); sol.C(grid.C.I)] - du;
        [sol.Phi(grid.Phi.I), sol.C(grid.C.I)] = ...
            split(u, grid.Phi.sz-2, grid.C.sz-2);
        sol = boundaries(sol, beta);
    end
    % Hessian matrix
    function H = hessian(sol)
        S11 = neumanns(grid.Phi, [1 0], [0 -1], [0 1]);
        S12 = dirichlets(grid.C, [-1 0], -1./sol.C(2, 2:end-1));
        
        S21 = dirichlets(grid.Phi, [-1 0], -exp(-sol.Phi(2, 2:end-1)));
        S22 = neumanns(grid.C, [0 -1], [0 1]);

        H = [L.Phi * S11, L.C * S12; ...
             L.Phi * S21, L.C * S22];
    end
    function S = dirichlets(grid, dir, vals)        
        K = find(~grid.I & shift(grid.I, dir));
        I = grid.I & ~shift(grid.I, -dir);
        J = find(I(2:end-1, 2:end-1));
        S = sparse(K, J, vals(:), grid.numel, nnz(grid.I));
    end
    function S = neumanns(grid, varargin)
        S = expand(grid.I);
        for k = 1:numel(varargin)
            d = varargin{k};
            S = neumann(grid, d) * S;
        end
    end
    step = @step_func;
end
