function main
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for the solver.
m = 6;
x = logspace(-1, 0, 1+2^m); 
y = 0; % logspace( -1, 0, 1+2^m); 

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
Ux = @(X, Y) cos(X);
Uy = @(X, Y) sin(Y);
C = @(X, Y) exp(X - Y);
L = @(X, Y) -exp(X - Y) .* (sin(X) - cos(X) - 0*cos(Y) + 0*sin(Y));

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, F, I] = laplacian(sz, X, Y, C(X, Y), L(X, Y));
if sz(1) > 1
    Bl = ~I & circshift(I, [-1 0]);
    Br = ~I & circshift(I, [+1 0]);
    [A, F] = boundary_dirichlet(A, F, Bl, X, Y, U);
    [A, F] = boundary_dirichlet(A, F, Br, X, Y, U);
end
if sz(2) > 1
    Bu = ~I & circshift(I, [0 -1]);
    Bd = ~I & circshift(I, [0 +1]); 
    [A, F] = boundary_dirichlet(A, F, Bu, X, Y, U);
    [A, F] = boundary_dirichlet(A, F, Bd, X, Y, U);
end
% [A, F] = boundary_neumann(A, F, ~I, Ux);

A = A(I, I); F = F(I);
fprintf('(%.3fs)\n', toc);
assert(nnz(A - A') == 0)

Vi = zeros(nnz(I), 1); % Initial guess.

fprintf('Construct iterative %s solver... ', type); tic;
[C, d] = jacobi(A, Vi, F); fprintf('(%.3fs)\n', toc);

fprintf('Apply %d iterations... ', iters); tic;
[Vf, residuals] = iterate(Vi, C, d, iters, type);  %#ok<NASGU>
fprintf('(%.3fs)\n', toc);

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
err = show(mat_file); 
fprintf('error: %e\n', err);
fprintf('max. eigenvalue: %f\n', max_eigs(C, 1));

function e = max_eigs(A, k)
e = eigs(A, k, 'lm', struct('disp', 0)); 
