% load
% v = []; for k = 1:numel(solutions), v(k) = solutions{k}.Vinf; end
% b = []; for k = 1:numel(solutions), b(k) = solutions{k}.beta; end
% loglog(b, v)
function solver
    tic;
    clc;
    filename = datestr(now, 'YYYYmmmdd_HHMMSS');
    sol0 = struct( ...
        'radius', logspace(0, 7, 400), ...
        'theta', linspace(0, pi, 50), ...
        'gamma', 0.01, ...
        'maxres', 1e-11, ...
        'iters', [2 5] ...    
    );
    betas = 10.^fliplr(-5:0.5:-1);
    sol0 = force(sol0); % initialize force solver
    sol0.Vinf = 1;
    solutions = cell(numel(betas), 1);
    alphas = [0, 0.5];
    
    for k = 1:numel(betas)
        sol = sol0;
        sol.beta = betas(k);
        for a = alphas
            sol.alpha = a;
            [sol, V, F] = steady(sol, sol.Vinf*(1+[0, 0.1]), 5);
            sol.Vinf = V(end);
        end
        solutions{k} = sol;
        pause(0);
%         pynotify('CEK', sprintf('%.1f%% completed.', 100*k/numel(betas)));
    end
    save(filename)
    show(filename)
end
