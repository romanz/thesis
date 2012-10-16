function plot_streamlines(sol)

s = streamfunc(sol);
[R, T] = ndgrid(s.grid.r, s.grid.t);
I = s.grid.r < 100;
R = R(I, :);
T = T(I, :);
Z = s.Psi(I, :);
X = R.*cos(T);
Y = R.*sin(T);
c = 0.001;
c = linspace(-c, c, 100);
contour(X, Y, Z, c);
axis equal