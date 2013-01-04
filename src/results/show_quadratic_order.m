function show_quadratic_order
load alpha=0.0_Du=0.5_zeta=6.0_Rmax=100.mat
Du = 0.5;
zeta = 6;
U1 = (Du*log(16) + zeta)/(1 + 2*Du);
U = repmat(U1*b, 1, size(Vn, 2));
E = abs(Vn - U);
clf;
set(gcf, 'DefaultAxesColorOrder', gray(5))
hold on;
plot(b, E, 'o', 'LineWidth', 2);
plot(b, b*0.17*(0.25.^[0:3]), ':', 'LineWidth', 2);

hold off;
set(gca, 'XScale', 'log', 'YScale', 'log')
legend({'32x32', '64x64', '128x128', '256x256'}, 'Location', 'SouthEast')
xlim([2e-2 0.1])
xlabel('Electric field magnitude (\beta)')
ylabel('Discretization error')
title('Quadratic convergence of steady-state velocity')
print -depsc2 quadratic_order
