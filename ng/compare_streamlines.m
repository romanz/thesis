%% Compute flows
clc;
clear;

gamma = exp(-1);
Vinf = 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));

beta = 1e-3;

% force_solver('res_stokes', 'w', beta, gamma, Vinf*beta, 1e4, ...
%     [120 40], 1000, 10*[0 1 0]);

opts = {'Rinf', 1e3, 'N', [60 15], 'cycles', 1000, 'iters', 2*[1 1 1], ...
    'filename', 'res_coupled', 'mode', 'w', 'quiet', true};

f = @(x) force_solver(beta, gamma, x*Vinf*beta, opts{:});
x = secant(f, [0.9 1.1], 3);
% f(1)


%%
load res_coupled;
f = 0;
figure(f+1); mesh(solVx);
figure(f+2); mesh(solVy);
figure(f+3); mesh(solP);
figure(f+4); mesh(solC);
figure(f+5); mesh(solPhi);

axs = [-1 1 0 1];

% figure(1); streamlines('res_stokes', 0, 2e-3, 20); axis equal; 
% axis(axs*7); title('Stokes Streamlines (numerical solution)')
% print -depsc2 stokes_numeric.eps
% 
% figure(2); 
% S = load('res_stokes', 'radius', 'theta', 'Vinf');
% gridPsi = init_grid(S.radius, S.theta);
% psi = @(r, theta) S.Vinf * sin(theta) .* (0.5*r - 0.75 + 0.25./(r.^2));
% contour(gridPsi.X .* cos(gridPsi.Y), gridPsi.X .* sin(gridPsi.Y), ...
%     psi(gridPsi.X, gridPsi.Y), linspace(0, 2e-3, 20));
% 
% axis equal; 
% axis(axs*7); title('Stokes Streamlines (theoretical solution)')
% print -depsc2 stokes_theory.eps

%
figure(f+6); 
streamlines('res_coupled', 0, beta, 100); axis equal; axis(3*axs); 
title('Coupled Flow Streamlines (numerical solution)')
print -depsc2 coupled_numeric.eps
%
figure(f+7); 
S = load('res_coupled', 'radius', 'theta', 'Vinf');
gridPsi = init_grid(S.radius, S.theta);
psi = @(r, theta) S.Vinf * sin(theta) .* (0.5*r - 0.75 + 0.25./(r.^2)) ...
                + 1.5*S.Vinf * sin(theta) .* (0.5 - 0.5./(r.^2));
contour(gridPsi.X .* cos(gridPsi.Y), gridPsi.X .* sin(gridPsi.Y), ...
    psi(gridPsi.X, gridPsi.Y), linspace(0, beta, 100));

axis equal; 
axis(axs*3); title('Coupled Flow Streamlines (theoretical solution)')
print -depsc2 coupled_theory.eps
%
figure(f+8)
Vslip = mean(solVy(1:2, :))';
plot(theta*180/pi, Vslip);
xlabel('\theta'); ylabel('V_{slip}')
