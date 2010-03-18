randn('state', 1);

iters = 10e3;
iter_type = 'Jacobi';
sz = 1+2.^([1 1]*5);
N = prod(sz);
x = linspace(-1, 1, sz(1));
y = linspace(-1, 1, sz(1));
[X, Y] = ndgrid(x, y);

[L, M, I] = laplacian(sz, X, Y);
L = dinv(M) * L;
% Rotational velocity field V(x, y)
H = [1;1]/2 * [0 1 0];
Vx =     average(Y, H);
Vy =    -average(X, H');
[VG] = gradient(sz, X, Y, Vx, Vy, 'centered');
A = L - VG;
U = @(X, Y) X.^2 + Y.^2;
f = 0 * X + 4;
f(~I) = NaN; % Not defined on the boundary

[A, f] = dirichlet(A, f, boundary(I, [-1 0]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [+1 0]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [0 -1]), X, Y, U);
[A, f] = dirichlet(A, f, boundary(I, [0 +1]), X, Y, U);
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(all(isnan(f(~I))) == 1);
assert(nnz(isnan(A(I, I))) == 0);
assert(nnz(isnan(f(I))) == 0);
A = A(I, I); 
f = f(I);

fprintf('Construct Jacobi iteration... '); tic;
[R, T, d] = jacobi(A, f); fprintf('(%.3fs)\n', toc);
fprintf('Maximal eigenvalue of T: ');
lambda = abs(max_eigs(T, 1));
fprintf('%.6f => %d iterations/decade\n', ...
    lambda, ceil(-1/log10(lambda)));

%% Iteration phase
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
