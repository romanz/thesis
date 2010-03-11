function main
addpath ../utils ../extrapolation
fprintf('\n');
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

%% Create grid for the solver.
m = 5;
x = logspace(-1, 0, 1+2*2^m); 
y = 1+linspace(-1, 0, 1+2^m); 

%% # of iterations
iters = 100e3;
iter_type = 'Jacobi';
% iter_type = 'RedBlack'; iters = iters / 2;
% iter_type = 'RRE'; iters = 4e3;
% iter_type = 'MPE'; iters = 4e3;

%% We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid(x, y);
sz = [numel(x) numel(y)];

%% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the diverence of C * grad(U).
% It is useful for solver's verification.
U = @(X, Y) sin(X) - cos(Y);
Ux = @(X, Y) cos(X);
Uy = @(X, Y) sin(Y);
C = @(X, Y) exp(X - Y);
L = @(X, Y) -exp(X - Y) .* (sin(X) - cos(X) - cos(Y) + sin(Y));

%% We actually solve the linear system: Av = f
fprintf('Compute Laplacian on %d x %d grid... ', numel(x), numel(y)); tic;
[A, f, I] = laplacian(sz, X, Y, C(X, Y), L(X, Y));

% Add boundary conditions for X:
if sz(1) > 1
    Bl = ~I & circshift(I, [-1 0]); % Left boundary (small X)
    Br = ~I & circshift(I, [+1 0]); % Right boundary (large X)
    
    [A, f] = boundary_dirichlet(A, f, Bl, X, Y, U);
    [A, f] = boundary_neumann(A, f, Br, [+1 0], X, Y, Ux, Uy);
%     [A, f] = boundary_dirichlet(A, f, Br, X, Y, U);
%     [A, f] = boundary_neumann(A, f, Bl, [-1 0], X, Y, Ux, Uy);
end
% Add boundary conditions for Y:
if sz(2) > 1
    Bd = ~I & circshift(I, [0 -1]); % Bottom boundary (small Y)
    Bu = ~I & circshift(I, [0 +1]); % Upper boundary (large Y)
    
%     [A, f] = boundary_dirichlet(A, f, Bu, X, Y, U);
    [A, f] = boundary_neumann(A, f, Bd, [0 -1], X, Y, Ux, Uy);
%     [A, f] = boundary_dirichlet(A, f, Bd, X, Y, U);
    [A, f] = boundary_neumann(A, f, Bu, [0 +1], X, Y, Ux, Uy);
end

%% Sanity checks for the linear system
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
% Verify that all entries are valid real numbers:
assert(nnz(isnan(A) | isinf(A)) == 0)
% Verify that the matrix is symmetric (since the original operator is self-adjoint):
assert(nnz(A - A') == 0)

%% Restrict the problem to interior variables and constuct iteration matrix
A = A(I, I); f = f(I);
fprintf('(%.3fs)\n', toc);

fprintf('Construct Jacobi iteration... '); tic;
[R, T, d] = jacobi(A, f); fprintf('(%.3fs)\n', toc);

fprintf('Maximal eigenvalue of T: ');
lambda = abs(max_eigs(T, 1));
fprintf('%.10f => %d iterations/decade\n', ...
    lambda, ceil(-1/log10(lambda)));

%% Iteration phase
randn('state', 1);
Vi = randn(nnz(I), 1); % Initial guess.

fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
if any(strcmpi(iter_type, {'MPE', 'RRE'}))
    cycle = 20; % NOTE: each iteration actually computes cycle+1 vectors.
    [Vf, residuals] = extrapolate(Vi, @(v) T*v + d, cycle, iters/cycle, iter_type);
else
    [Vf, residuals] = iterate(Vi, A, f, R, iters, iter_type); 
end
fprintf('(%.3fs)\n', toc);

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
err = show(mat_file); 
fprintf('error: %e\n', err);
% full(A), f
