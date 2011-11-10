% load
% v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
% b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
% loglog(b, v)
function solver
    tic;
    clc;
    sol = struct( ...
        'radius', logspace(0, 7, 400), ...
        'theta', linspace(0, pi, 50), ...
        'gamma', 1-.001, ...
        'alpha', 0, ...
        'maxres', 1e-11, ...
        'iters', [2 6] ...    
    );
    betas = 10.^(-4:1:0);
    sol = force(sol); % initialize force solver
    solutions = cell(numel(betas), 1);
    for k = 1:numel(betas)
        sol.beta = betas(k);
        sol = analytic(sol); % approximation
        alphas = [0]; % XXX
        for a = alphas
            sol.alpha = a;
            [sol, V, F] = steady(sol, sol.Vinf + [0, 0.001], 10);
        end
        solutions{k} = sol;
        pause(0);
    end
    save
    tone(300, 0.1)
    b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
    v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
    plot(log10(b), log10([v; 0.042*b.^3 + 5e-4*b]), '.-')
end

