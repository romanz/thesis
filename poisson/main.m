function main
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for the solver.
m = 6;
x = logspace(-1, 0, 1+2^m); 
y = logspace( -1, 0, 1+2^m); 

% # of iterations
iters = 25e3;
% type = 'RedBlack';
type = 'Jacobi';

% We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid(x, y);
sz = [numel(x) numel(y)];
% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It is useful for solver's verification.
U = @(X, Y) sin(X) - cos(Y);
C = @(X, Y) exp(X - Y);
L = @(X, Y) -exp(X - Y) .* (sin(X) - cos(X) - cos(Y) + sin(Y));

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, F, I] = laplacian(sz, X, Y, C(X, Y), L(X, Y));
[A, F] = boundary_dirichlet(A, F, ~I, X, Y, U);

A = A(I, I); F = F(I);
fprintf('(%.3fs)\n', toc);

V0 = U(X, Y); % The ideal continuous solution.
Vi = zeros(nnz(I), 1); % Initial guess zeros(size(V0(I))

fprintf('Construct iterative %s solver... ', type); tic;
[C, d] = jacobi(A, Vi, F); fprintf('(%.3fs)\n', toc);

fprintf('Apply %d iterations... ', iters); tic;
[Vf, residuals] = iterate(Vi, C, d, iters, type); 
fprintf('(%.3fs)\n', toc);

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
show(mat_file)
assert(nnz(A - A') == 0)
