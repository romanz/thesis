clear;
clf;
load ~/MATs/sol_beta=1.000e-01_[129x129]_Rmax=100.0_Du=1.00_zeta=10.00_alpha=0.50.mat
% op = sol.C;
% Z = regrid(op);
% g = op.grid;
% h = [1 1]/2;
% R = convn(g.R, h, 'valid');
% T = convn(g.T, h, 'valid');
% z = convn(Z, h, 'valid');
% x = R .* cos(T);
% y = R .* sin(T);
% hold on
% c = 1 + linspace(-1, 1, 100);
% contourf(x, y, z, c, '-');
% contourf(x, -y, z, c, '-');
% t = linspace(0, 2*pi, 100);
% fill(cos(t), sin(t), [1 1 1]*0.5)
% axis([-1 1 -1 1]*4)
% axis equal
% colorbar
% return 

% s = streamfunc(sol);
% I = find(s.grid.Psi.r < 10);
% R = s.grid.Psi.R(I, :);
% T = s.grid.Psi.T(I, :);
% z = s.Psi(I, :);
% x = R .* cos(T);
% y = R .* sin(T);
% c = linspace(-1, 1, 61)*3 * sol.beta;
% hold on;
% contour(x, y, z, c)
% contour(x, -y, z, c)
% t = linspace(0, 2*pi, 100);
% fill(cos(t), sin(t), [1 1 1]*0.5)
% axis([-1 1 -1 1]*3)
% axis equal
% return

clf;
hold on;
op = sol.Phi;
g = op.grid;
I = g.r < 6;
R = g.R(I, :);
T = g.T(I, :);
X = R .* cos(T);
Y = R .* sin(T);
Z = regrid(op);
Z = Z(I, :);
Er = diff(Z(1:2, :), 1, 1) ./ diff(g.r(1:2));
Er = (Er(1:end-1) + Er(2:end)) / 2;
Et = diff(mean(Z(1:2, :)), 1, 2) ./ diff(g.t(1:2));
whos Er Et
x = convn(X, [1 1; 1 1]/4, 'valid');
y = convn(Y, [1 1; 1 1]/4, 'valid');
z = convn(Z, [1 1; 1 1]/4, 'valid');
z = z - z(1, end);
c = (-1:0.01:1)*sol.beta*10;
contourf(x, y, z, c)
contourf(x, -y, z, c)
t = linspace(0, 2*pi, 1000);
fill(cos(t), sin(t), [1 1 1]*0.5)
quiver(cos(sol.grid.t), sin(sol.grid.t), Er(:), Et(:), 'k')
quiver(cos(sol.grid.t), -sin(sol.grid.t), Er(:), -Et(:), 'k')
axis([-1 1 -1 1]*3)
axis equal

