function streamlines_test
clc
load results
% [X, Y] = ndgrid(logspace(0, 1), linspace(0, pi));
% psi = Vinf * sin(Y)/2 .* (X - 3./(2*X) + 1./(2*X.^2));
% Xp = X .* cos(Y);
% Yp = X .* sin(Y);
% clf;
% contour(Xp, Yp, psi, linspace(0, .05, 100))
% axis equal
% axis([-6 6 0 6])

% radial velocity component
gridPsi = init_grid(radius, theta);
% D1 = grad(gridPsi.sz, 2);
% D1 = spdiag(1./(D1 * gridPsi.Y(:))) * D1 * spdiag(sin(gridPsi.Y));
% D1 = spdiag(1./(gridVx.X(:, 2:end-1) .* sin(gridVx.Y(:, 2:end-1)))) * D1;
% 
% D2 = grad(gridPsi.sz, 1);
% D2 = spdiag(1./(D2 * gridPsi.X(:))) * D2 * spdiag(sin(gridPsi.X));
% D2 = -spdiag(1./(gridVy.X(2:end-1, :))) * D2;

D1 = grad(gridVx.X(:, 2:end-1), 1);
D1 = D1 * spdiag(gridVx.X(:, 2:end-1) .^ 2); 
D1 = spdiag(1./(gridC.X(gridC.I) .^ 2)) * D1;
d1 = D1 * col(solVx(:, 2:end-1));

D2 = grad(gridVy.Y(2:end-1, :), 2);
D2 = D2 * spdiag(sin(gridVy.Y(2:end-1, :)));
D2 = spdiag(1./(sin(gridC.Y(gridC.I)) .* gridC.X(gridC.I))) * D2;
d2 = D2 * col(solVy(2:end-1, :));

div = d1 + d2;
norm(div(:), inf)

Vx = solVx .* sin(gridVx.Y) .* gridVx.X.^2;
Vx = Vx(:, 2:end-1);
Vy = solVy .* sin(gridVy.Y) .* gridVy.X;
Vy = Vy(2:end-1, :);
V = [Vx(:); Vy(:)];

d = [grad(gridVx.X(:, 2:end-1), 1) grad(gridVy.Y(2:end-1, :), 2)] * V;
norm(d, inf)

G1 = +grad(gridPsi.Y, 2);
G2 = -grad(gridPsi.X, 1);
G = [G1; G2];
I = gridPsi.I;
I = I | shift(I, [1 0]);
P = expand(I);
Psi = P * ((G*P) \ V);
Psi = Psi - Psi(1);
norm(G * Psi - V, inf)
Psi = reshape(Psi, gridPsi.sz);
W = (gridPsi.X .* sin(gridPsi.Y));
Psi(~~W) = -Psi(~~W) ./ W(~~W);

[~, ~, Psi] = stokes_solution(Vinf, gridPsi.X, gridPsi.Y);
contour(gridPsi.X .* cos(gridPsi.Y), gridPsi.X .* sin(gridPsi.Y), ...
    Psi, linspace(0, max(Psi(:)), 500));
R = 3;
axis equal
axis([-R R 0 R]*2)
% D = [D1; D2];
% V = [col(solVx(:, 2:end-1)); col(solVy(2:end-1, :))];
% 
% D * V
% 
% psi = reshape(psi, gridPsi.sz);
% % mesh(theta, radius, psi)
% 
% Xpsi = gridPsi.X .* cos(gridPsi.Y);
% Ypsi = gridPsi.X .* sin(gridPsi.Y);
% 
% psi0 = Vinf * sin(gridPsi.Y)/2 .* (gridPsi.X - 3./2 + 1./(2*gridPsi.X.^2));
% 
% clf
% I = gridPsi.x < 1.5;
% mesh(gridPsi.X(I, :), gridPsi.Y(I, :), psi(I, :))
% 
% % contour(Xpsi, Ypsi, psi, linspace(0, .1, 200)); axis equal; axis([-1 1 0 2] * 10)
% 
% return

function D = grad(X, dim)
sz = size(X);
I = true(sz);
dir = (1:numel(sz)) == dim;

Kn = find(shift(I, -dir));
Kp = find(shift(I, +dir));
D = sparse(1:numel(Kp), Kp, 1, numel(Kp), numel(I)) - ...
    sparse(1:numel(Kn), Kn, 1, numel(Kn), numel(I));

D = spdiag(1 ./ (D * X(:))) * D;