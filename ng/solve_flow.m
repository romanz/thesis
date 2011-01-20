%% Compute flows
clear;
fprintf('\n');

gamma = exp(-1);
Vinf = 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));

beta = .1

% force_solver('res_stokes', 'w', beta, gamma, Vinf*beta, 1e4, ...
%     [120 40], 1000, 10*[0 1 0]);

opts = {'Rinf', 1e3, 'N', [60 15], 'cycles', 100, 'iters', 2*[1 1 1], ...
    'filename', 'res_coupled', ...
    'mode', 'rw', 'quiet', true};

f = @(x) force_solver(beta, gamma, x*Vinf*beta, opts{:});
x = secant(f, [0.95 1.05], 1); % Secant method for zero-crossing
% f(0.886280)
% secant(f, [0.61573 0.6158], 100); % Secant method for zero-crossing
