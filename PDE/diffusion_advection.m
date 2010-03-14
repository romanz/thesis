iters = 10e3;
iter_type = 'RedBlack';
sz = 1+2.^[3 3];
N = prod(sz);
x = linspace(-1, 1, sz(1));
y = linspace(-1, 1, sz(1));
[X, Y] = ndgrid(x, y);

[L, M, I] = laplacian(sz, X, Y);
L = dinv(M) * L;

% Rotational velocity field V(x, y)
Vx =    Y;
Vy = -3*X;
[Gx, Gy] = gradient(sz, X, Y);
A = L - (spdiag(Vx)*Gx + spdiag(Vy)*Gy);
U = @(X, Y) X.^2 + 2*Y.^2;
f = 10 * X .* Y + 6;

[A, f] = dirichlet(A, f, boundary(I, [-1 0]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [+1 0]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [0 -1]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [0 +1]), X, Y, U);
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
% assert(nnz(f(~I)) == 0); % XXX
A = A(I, I); 
f = f(I);

fprintf('Construct Jacobi iteration... '); tic;
[R, T, d] = jacobi(A, f); fprintf('(%.3fs)\n', toc);

%% Iteration phase
randn('state', 1);
Ui = randn(nnz(I), 1); % Initial guess.

fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
[Uf, residuals] = iterate(Ui, A, f, R, iters, iter_type); 
fprintf('(%.3fs)\n', toc);

%% Save and show the results.
figure(1); clf;
mat_file = 'results.mat';
save(mat_file)
err = show(mat_file); 
fprintf('error: %e\n', err);

figure(2); clf;
quiver(X, Y, Vx, Vy);
