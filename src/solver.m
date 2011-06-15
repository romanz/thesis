function sol = solver
    clc;
    tic;
    sol = struct( ...
        'radius', logspace(0, 3, 60), ...
        'theta', linspace(0, pi, 40), ...
        'gamma', 10, ...
        'alpha', 0, ...
        'maxres', 1e-9, ...
        'iters', [5, 20] ...
    );
    sol = force(sol); % initialize force solver

    betas = 1; % logspace(0, 1, 20);
    sol.Vinf = 1;

    alphas = 0; [0 0.1];
    for a = alphas
        sol.alpha = a;
        [sol, V, F] = steady(sol, betas(1), sol.Vinf*[1, 0.9], 200);
    end
    save
end

function [sol, V, F] = steady(sol, betas, v, iters)

    logger('steady', 'alpha = %f', sol.alpha)
    function f = force_func(b, u) % Total force for specified (beta, Vinf)
        [sol] = force(sol, b, u);
        f = sol.force.total;
    end

    force_func(betas, v(1));
    % Last beta, for steady-state velocity secant method.
    b = betas(end);
    logger('steady', 'secant (beta = %e)', b)
    secant_step = secant(@(u) force_func(b, u), v);
    for k = 1:iters
        [V(k), F(k)] = secant_step();
        log('steady', '(%d/%d) V = %e \t F = %e', k, iters, V(k), F(k));
        plot([V; F]'); drawnow;
    end
    
end
