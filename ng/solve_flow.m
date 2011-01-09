%% Compute flows
clear;
fprintf('\n');

gamma = exp(-1);
Vinf = 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));

beta = .2;

% force_solver('res_stokes', 'w', beta, gamma, Vinf*beta, 1e4, ...
%     [120 40], 1000, 10*[0 1 0]);

opts = {'Rinf', 1e3, 'N', [60 15], 'cycles', 1000, 'iters', 2*[1 1 1], ...
    'filename', 'res_coupled', ...
    'mode', '', 'quiet', true};

f = @(x) force_solver(beta, gamma, x*Vinf*beta, opts{:});
x = secant(f, [0.9 1.1], 3); % Secant method for zero-crossing
% f(1)
