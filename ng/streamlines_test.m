function streamlines_test
clc
load results
[X, Y] = ndgrid(logspace(0, 1), linspace(0, pi));
psi = Vinf * sin(Y)/2 .* (X - 3./(2*X) + 1./(2*X.^2));
Xp = X .* cos(Y);
Yp = X .* sin(Y);
clf
contour(Xp, Yp, psi, linspace(0, .05, 50))
axis equal