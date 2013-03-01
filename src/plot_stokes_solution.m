clear all;
clf;

r = linspace(1, 10, 513);
t = linspace(0, pi, 513);
[R, T] = ndgrid(r, t);
beta = 0.1;

Du = 1;
zeta = 10;
U = beta * (Du*log(16) + zeta)/(1 + 2*Du);
Psi = U * 0.5 * (R.^(-2) - 0*R) .* sin(T);

z = Psi;
z(:, [1 end]) = 1e-10;
x = R .* cos(T);
y = R .* sin(T);
c = linspace(-1, 1, 101)*beta*10;
hold on;
contour(x, y, z, c)
contour(x, -y, z, c)
t = linspace(0, 2*pi, 100);
fill(cos(t), sin(t), [1 1 1]*0.5)
axis([-1 1 -1 1]*2.5)
axis equal
title(sprintf('Streamlines (\\Psi): \\beta = %.2f', beta))
print('-depsc2', 'LinearStokes_Psi-uniformflow.eps')
