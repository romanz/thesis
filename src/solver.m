function sol = solver
    tic;
    sol = struct( ...
        'radius', logspace(0, 3, 60), ...
        'theta', linspace(0, pi, 40), ...
        'gamma', 10, ...
        'alpha', 0, ...
        'maxres', 1e-9, ...
        'iters', [1, 20] ...
    );
    sol = force(sol); % initialize force solver

    betas = 1; % logspace(0, 1, 20);
    
    sol = steady(sol, betas(1), [1, 0.9], 5);
    sol.alpha = 0.5;
    sol = steady(sol, betas(1), sol.Vinf*[1, 0.9], 5);
end

function [sol] = steady(sol, betas, v, iters)

    log('steady', 'alpha = %f', sol.alpha)
    function f = force_func(b, u) % Total force for specified (beta, Vinf)
        [sol] = force(sol, b, u);
        f = sol.force.total;
    end

    force_func(betas, v(1));
    % Last beta, for steady-state velocity secant method.
    b = betas(end);
    log('steady', 'secant (beta = %e)', b)
    secant_step = secant(@(u) force_func(b, u), v);
    for k = 1:iters
        log('steady', 'secant %d/%d', k, iters)
        secant_step();
    end
    
end
