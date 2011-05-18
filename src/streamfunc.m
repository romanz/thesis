function sol = streamfunc(sol)

grid = sol.grid;
grid.Psi = init_grid(grid.radius, grid.theta);
% d(r^2 sint Vr)/dr + d(r sint Vt)/dt = 0

% d(Vx)/dx + d(Vy)/dy = 0
% Vx = r^2 sint Vx = dPsi/dy
% Vy = r sint Vy = - dPsi/dx

Vx = sol.Vx .* sin(grid.Vx.Y) .* grid.Vx.X.^2;
Vx = Vx(:, 2:end-1);
Vy = sol.Vy .* sin(grid.Vy.Y) .* grid.Vy.X;
Vy = Vy(2:end-1, :);
V = [Vx(:); Vy(:)];

div = [grad(grid.Vx.X(:, 2:end-1), 1) grad(grid.Vy.Y(2:end-1, :), 2)] * V;
assert(norm(div) / norm(V) < 1e-10);

G1 = +grad(grid.Psi.Y, 2);
G2 = -grad(grid.Psi.X, 1);
G = [G1; G2]; % [d/dy; -d/dx] * Psi = V = [Vx; Vy]

I = grid.Psi.I | shift(grid.Psi.I, [1 0]); % I = 0 iff Psi = 0
Q = expand(I); % G*(Q*Psi) = V

Psi = Q * ((G*Q) \ V);
assert(norm(G * Psi - V) / norm(V) < 1e-10);
Psi = reshape(Psi, grid.Psi.sz);

w = (grid.Psi.X .* sin(grid.Psi.Y));
Psi(I) = Psi(I) ./ w(I);

sol.Psi = Psi;
sol.grid = grid;

function D = grad(X, dim)
sz = size(X);
I = true(sz);
dir = (1:numel(sz)) == dim;

Kn = find(shift(I, -dir));
Kp = find(shift(I, +dir));
D = sparse(1:numel(Kp), Kp, 1, numel(Kp), numel(I)) - ...
    sparse(1:numel(Kn), Kn, 1, numel(Kn), numel(I));

D = spdiag(1 ./ (D * X(:))) * D;