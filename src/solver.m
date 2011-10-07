function solver
    tic;
    clc;
    sol = struct( ...
        'radius', logspace(0, 3, 75), ...
        'theta', linspace(0, pi, 25), ...
        'gamma', 2, ...
        'alpha', 0, ...
        'maxres', 1e-12, ...
        'iters', [2 10] ...    
    );
    betas = 10.^(-2:0.01:0);
    sol = force(sol); % initialize force solver
    solutions = cell(numel(betas), 1);
    for k = 1:numel(betas)
        sol.beta = betas(k);
        sol = analytic(sol); % approximation
        alphas = [0 0.5]; % XXX
        for a = alphas
            sol.alpha = a;
            [sol, V, F] = steady(sol, sol.Vinf*[1, 0.9], 5);
        end
        solutions{k} = sol;
        pause(0);
    end
    save
    tone(300, 0.1)
end
