% load
% v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
% b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
% loglog(b, v)
function solver
    tic;
    clc;
    sol = struct( ...
        'radius', logspace(0, 4, 95), ...
        'theta', linspace(0, pi, 25), ...
        'gamma', 1.001, ...
        'alpha', 0, ...
        'maxres', 1e-12, ...
        'iters', [2 10] ...    
    );
    betas = 10.^(-3:0.1:0);
    sol = force(sol); % initialize force solver
    solutions = cell(numel(betas), 1);
    for k = 1:numel(betas)
        sol.beta = betas(k);
        sol = analytic(sol); % approximation
        alphas = [0]; % XXX
        for a = alphas
            sol.alpha = a;
            [sol, V, F] = steady(sol, sol.Vinf + [0, 0.001], 5);
        end
        solutions{k} = sol;
        pause(0);
    end
    save
    tone(300, 0.1)
end

