function streamlines_test
clc
load results
[X, Y] = ndgrid(logspace(0, 1), linspace(0, pi));
psi = Vinf * sin(Y)/2 .* (X - 3./(2*X) + 1./(2*X.^2));
Xp = X .* cos(Y);
Yp = X .* sin(Y);
clf;
contour(Xp, Yp, psi, linspace(0, .05, 100))
axis equal
axis([-6 6 0 6])

% radial velocity component
Dx = spdiag(1 ./ (gridVx.X.^2 .* sin(gridVx.Y))) * grad(gridVx.sz, 2);
Dy = spdiag(1 ./ (gridVx.X.^2 .* sin(gridVx.Y))) * grad(gridVx.sz, 2);

function D = grad(sz, dim)
I = true(sz);
dir = (1:numel(sz)) == dim;

Kn = find(shift(I, -dir));
Kp = find(shift(I, +dir));
D = sparse(1:numel(Kp), Kp, 1, numel(Kp), numel(I)) - ...
    sparse(1:numel(Kn), Kn, 1, numel(Kn), numel(I));
