
sz = 1+2.^[3 3];
N = prod(sz);
x = linspace(-1, 1, sz(1));
y = linspace(-1, 1, sz(1));
[X, Y] = ndgrid(x, y);

[L, M, I] = laplacian(sz, X, Y);
L = dinv(M) * L;

% Rotational velocity field: V(x, y) = (-y, x)
Vx = spdiag( Y);
Vy = spdiag(-X);
[Gx, Gy] = gradient(sz, X, Y);
A = L - (Vx*Gx + Vy*Gy);
U = @(X, Y) X.^2 + Y.^2;
f = 0*X + 0*Y + 4;

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

iters = 10e3;
iter_type = 'Jacobi';
fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
[Uf, residuals] = iterate(Ui, A, f, R, iters, iter_type); 
fprintf('(%.3fs)\n', toc);

u = A \ f;
U1 = U(X, Y);
U1(I) = u;
U2 = U(X, Y);
U2(I) = Uf;
U2 - U1
