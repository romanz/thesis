function test_laplacian

opts = {'Rinf', 1e3, 'N', [200 100], 'cycles', 100, 'iters', 2*[1 1 1], ...
    'filename', 'res_coupled', ...
    'mode', 'rw', 'quiet', true};

opts = struct(opts{:});
radius = logspace(0, log10(opts.Rinf), opts.N(1)).';
theta = linspace(0, pi, opts.N(2)).';
[center, interior, xstag, ystag] = grids(radius, theta);

beta = 1;
[A, a, A1] = maxwell(center, beta);
[B, b, B1] = advection(center);
whos A B a b A1 B1

rhs = [a; b];
M = [A A1; B1 B];
tic;
u = zeros(size(rhs));
for i = 1:10
    r = rhs - M*u;
    du = M \ r;
    u = u + du;
    disp(norm(r(:), inf))
end
toc;
whos u

u1 = u(1:end/2);
u2 = u((1+end/2):end);
u1 = reshape(u1, center.sz-2);
u2 = reshape(u2, center.sz-2);
X = center.X(2:end-1, 2:end-1);
Y = center.Y(2:end-1, 2:end-1);
plot([mean([u1(1, :); u2(1, :)])', beta*0.75*cos(center.y(2:end-1))])
% figure(1); mesh(u1 - (cos(Y).*(X - 0.25./X.^2))) % Phi
% figure(2); mesh(u2 - (cos(Y).*(0.75./X.^2))) % C
% u = P*u + q;
% u = reshape(u, center.sz-2);
% figure(1); mesh(u); 
% u0 = cos(center.Y(2:end-1, 2:end-1)) .* (center.X(2:end-1, 2:end-1) );
% figure(2); mesh(u0); 
% err = u - u0;
% err = err(1:10, :);
% err = err(:);
% disp(norm(err, 2) / sqrt(numel(u)))
% disp(norm(err, inf))
end

function [A, rhs, A1] = maxwell(grid, b)
L = laplacian(grid.I, grid.X, grid.Y);
S = neumann(grid, [0 -1], 1) * neumann(grid, [0 +1], 1);
[Pe, Qe] = neumann(grid, [1 0]);
P = S * Pe * expand(grid.I);
q = S * (Qe * (b * cos(grid.y(2:end-1))) );

Q1 = dirichlet(grid, [-1 0]);
J = false(grid.sz-2);
J(1, :) = 1;
J1 = expand(J);
A1 = L * S * Q1 * J1';

rhs = -L * q;
A = L * P;
end

function [A, rhs, A1] = advection(grid)
L = laplacian(grid.I, grid.X, grid.Y);
S = neumann(grid, [0 -1], 1) * neumann(grid, [0 +1], 1);
[Q] = dirichlet(grid, [1 0]);
P = S * expand(grid.I);
q = S * (Q * zeros(size(grid.y(2:end-1))) );

Q1 = dirichlet(grid, [-1 0]);

J = false(grid.sz-2);
J(1, :) = 1;
J1 = expand(J);
A1 = L * S * Q1 * J1';
rhs = -L * q;
A = L * P;
end
