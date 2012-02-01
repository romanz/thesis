% load
% v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
% b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
% loglog(b, v)
function solver
    tic;
    clc;
    filename = datestr(now, 'YYYYmmmdd-HHMMSS');
    sol = struct( ...
        'radius', logspace(0, 7, 400), ...
        'theta', linspace(0, pi, 50), ...
        'gamma', 1, ...
        'maxres', 1e-11, ...
        'iters', [2 5] ...    
    );
    betas = 10.^([-2:0.1:1]);
    sol = force(sol); % initialize force solver
    sol.Vinf = 0.1;
    solutions = cell(numel(betas), 1);
    alphas = [0];
    for k = 1:numel(betas)
        sol.beta = betas(k);
        for a = alphas
            sol.alpha = a;
            [sol, V, F] = steady(sol, sol.Vinf*(1+[-0.1, 0.1]), 5);
            sol.Vinf = V(end);
        end
        solutions{k} = sol;
        pause(0);
        pynotify('CEK', sprintf('%.1f%% completed.', 100*k/numel(betas)));
    end
    save(filename)
    show(filename)
end
