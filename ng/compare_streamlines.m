%% Compute flows
clear;

gamma = exp(-1);
Vinf = 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));

beta = 1e-3;
force_solver('res_stokes', 'w', beta, gamma, Vinf*beta, 1e4, ...
    [120 40], 1000, 10*[0 1 0]);

force_solver('res_coupled', 'w', beta, gamma, Vinf*beta, 1e4, ...
    [120 40], 1000, 10*[1 1 1]);

%% Print figures
axs = [-1 1 0 1]*7;

figure(1); streamlines('res_stokes', 2e-3, 20); axis equal; 
axis(axs); title('Stokes Streamlines (numerical solution)')
print -depsc2 stokes_numeric.eps

figure(2); 
S = load('res_stokes', 'radius', 'theta', 'Vinf');
gridPsi = init_grid(S.radius, S.theta);
psi = @(r, theta) S.Vinf * sin(theta) .* (0.5*r - 0.75 + 0.25./(r.^2));
contour(gridPsi.X .* cos(gridPsi.Y), gridPsi.X .* sin(gridPsi.Y), ...
    psi(gridPsi.X, gridPsi.Y), linspace(0, 2e-3, 20));

axis equal; 
axis(axs); title('Stokes Streamlines (theoretical solution)')
print -depsc2 stokes_theory.eps

figure(3); streamlines('res_coupled', 2e-3, 20); axis equal; 
axis(axs); title('Coupled Flow Streamlines (numerical solution)')
print -depsc2 coupled_numeric.eps

figure(4); 
S = load('res_coupled', 'radius', 'theta', 'Vinf');
gridPsi = init_grid(S.radius, S.theta);
psi = @(r, theta) S.Vinf * sin(theta) .* (0.5*r - 0.75 + 0.25./(r.^2)) ...
                + 1.5*S.Vinf * sin(theta) .* (0.5 - 0.5./(r.^2));
contour(gridPsi.X .* cos(gridPsi.Y), gridPsi.X .* sin(gridPsi.Y), ...
    psi(gridPsi.X, gridPsi.Y), linspace(0, 2e-3, 20));

axis equal; 
axis(axs); title('Coupled Flow Streamlines (theoretical solution)')
print -depsc2 coupled_theory.eps
