%% Compute flows
clear;

gamma = exp(-1);
Vinf = 2*log((gamma^0.25 + gamma^-0.25) / (2*gamma^0.25));
R = 3;

beta = 1e-3;
force_solver('res_stokes', 'w', beta, gamma, Vinf*beta, 1e4, ...
    [120 40], 1000, 10*[0 1 0]);

force_solver('res_coupled', 'w', beta, gamma, Vinf*beta, 1e4, ...
    [120 40], 1000, 10*[1 1 1]);

%% Print figures
figure(1); streamlines('res_stokes', 0, 2e-3, 20); axis equal; 
axis([-R R 0 R]*2); title('Stokes Streamlines (numerical solution)')
print -depsc2 stokes_numeric.eps

figure(2); streamlines('res_stokes', 1, 2e-3, 20); axis equal; 
axis([-R R 0 R]*2); title('Stokes Streamlines (theoretical solution)')
print -depsc2 stokes_theory.eps

figure(3); streamlines('res_coupled', 0, 2e-3, 20); axis equal; 
axis([-R R 0 R]*2); title('Coupled Flow Streamlines (numerical solution)')
print -depsc2 coupled_numeric.eps
