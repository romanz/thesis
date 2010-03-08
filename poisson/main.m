function main
fprintf('\n');
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for the solver.
m = 7;
x = logspace(-1, 0, 1+2^m); 
y = logspace(-1, 0, 1+2^m); 

% # of iterations
iters = 80e3;
% iter_type = 'RedBlack';
iter_type = 'Jacobi';

% We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid(x, y);
sz = [numel(x) numel(y)];
% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the diverence of C * grad(U).
% It is useful for solver's verification.
U = @(X, Y) sin(X) - cos(Y);
Ux = @(X, Y) cos(X);
Uy = @(X, Y) sin(Y);
C = @(X, Y) exp(X - Y);
L = @(X, Y) -exp(X - Y) .* (sin(X) - cos(X) - cos(Y) + sin(Y));

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, f, I] = laplacian(sz, X, Y, C(X, Y), L(X, Y));
if sz(1) > 1
    Bl = ~I & circshift(I, [-1 0]);
    Br = ~I & circshift(I, [+1 0]);
    [A, f] = boundary_dirichlet(A, f, Bl, X, Y, U);
%     [A, f] = boundary_neumann(A, f, Br, [+1 0], X, Y, Ux, Uy);
    [A, f] = boundary_dirichlet(A, f, Br, X, Y, U);
%     [A, f] = boundary_neumann(A, f, Bl, [-1 0], X, Y, Ux, Uy);
end
if sz(2) > 1
    Bd = ~I & circshift(I, [0 -1]); 
    Bu = ~I & circshift(I, [0 +1]);
    [A, f] = boundary_dirichlet(A, f, Bu, X, Y, U);
%     [A, f] = boundary_neumann(A, f, Bd, [0 -1], X, Y, Ux, Uy);
    [A, f] = boundary_dirichlet(A, f, Bd, X, Y, U);
%     [A, f] = boundary_neumann(A, f, Bu, [0 +1], X, Y, Ux, Uy);
end
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
% All entries are valid real numbers
assert(nnz(isnan(A) | isinf(A)) == 0)
% The matrix should be symmetric (since the original operator is self-adjoint)
assert(nnz(A - A') == 0)
A = A(I, I); f = f(I);
fprintf('(%.3fs)\n', toc);

randn('state', 1);
Vi = randn(nnz(I), 1); % Initial guess.

fprintf('Construct iterative %s solver... ', iter_type); tic;
[R, T, d] = jacobi(A, Vi, f); fprintf('(%.3fs)\n', toc);

fprintf('Apply %d iterations... ', iters); tic;
% [Vf, residuals] = iterate(Vi, A, f, R, iters, iter_type);  %#ok<NASGU>
[Vf, residuals] = extrapolate(Vi, @(v) T*v + d, 20, 200, 'RRE');
% [Vf, residuals] = extrapolate(Vi, @(v) T*v + d, 20, 200, 'MPE');
fprintf('(%.3fs)\n', toc);

lambda = abs(max_eigs(T, 1));
fprintf('Maximal eigenvalue of T: %.10f\n', lambda);

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
err = show(mat_file); 
fprintf('error: %e\n', err);
% full(A), f

function e = max_eigs(A, k)
if nnz(A) > 0
    e = eigs(A, k, 'lm', struct('disp', 0)); 
else
    e = 0;
end
