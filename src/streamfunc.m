function sol = streamfunc(sol)

grid = sol.grid;
grid.Psi = init_grid(grid.r, grid.t);
% d(r^2 sint Vr)/dr + d(r sint Vt)/dt = 0

% d(Vr)/dx + d(Vt)/dy = 0
% Vr = r^2 sint Vr = dPsi/dy
% Vt = r sint Vt = - dPsi/dx

Vr = regrid(sol.Vr);
Vt = regrid(sol.Vt);

Vr = Vr .* sin(grid.Vr.T) .* grid.Vr.R.^2;
Vr = Vr(:, 2:end-1);
Vt = Vt .* sin(grid.Vt.T) .* grid.Vt.R;
Vt = Vt(2:end-1, :);
V = [Vr(:); Vt(:)];

div = [grad(grid.Vr.R(:, 2:end-1), 1) grad(grid.Vt.T(2:end-1, :), 2)] * V;
assert(norm(div) / norm(V) < 1e-10);

G1 = +grad(grid.Psi.T, 2);
G2 = -grad(grid.Psi.R, 1);
G = [G1; G2]; % [d/dy; -d/dx] * Psi = V = [Vr; Vt]

I = grid.Psi.I | shift(grid.Psi.I, [1 0]); % I = 0 iff Psi = 0
Q = expand(I); % G*(Q*Psi) = V

Psi = Q * ((G*Q) \ V);
assert(norm(G * Psi - V) / norm(V) < 1e-10);
Psi = reshape(Psi, grid.Psi.sz);

w = (grid.Psi.R .* sin(grid.Psi.T));
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