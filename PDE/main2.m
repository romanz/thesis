fprintf('\n');
%% Simulate and solve Diffusion-Advection Equation on a 2D grid.
% The equation to solve is: D Laplacian(u) - V*Gradient(u) = f.
% Central difference is used for Laplacian. 
% Upstream difference is used for Gradient.
% Jacobi iteration is used for relaxation.
randn('state', 1);
iters = 1e3;
iter_type = 'MPE';
sz = 1+2.^([1 1]*7);

%% Create the grid
N = prod(sz);
x = linspace(-1, 1, sz(1));
y = linspace(-1, 1, sz(1));
[X, Y] = ndgrid(x, y);

%% L is the discrete Laplacian operator
[L, M, I] = laplacian(sz, X, Y);
L = dinv(M) * L;

%% Velocity field (Vx, Vy) definition on a staggered grid
H = [1;1]/2 * [0 1 0]; % Average on X, remove Y's boundary
Vx =     average(Y, H ); % Evaluated at left and right edges.
Vy =    -average(X, H'); % Evaluated at top and bottom edges.
%% VG is the discrete advection operator
[VG] = gradient(sz, X, Y, Vx, Vy, 'upstream');

%% The operator combined operator
A = L - VG; 
U = @(X, Y) X.^2 - Y.^2;
f = -4*X .* Y + 0; % Right-Hand Side
f(~I) = NaN; % Not defined on the boundary

%% Substitute Boundary Conditions
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
% Restrict the operator and RHS on the interior variables
A = A(I, I); 
f = f(I);

%% Apply relaxation to solve the problem
fprintf('Construct Jacobi iteration... '); tic;
[R, T, d] = jacobi(A, f); fprintf('(%.3fs)\n', toc);
fprintf('Maximal eigenvalue of T: ');
lambda = abs(max_eigs(T, 1));
fprintf('%.6f => %d iterations/decade\n', ...
    lambda, ceil(-1/log10(lambda)));

%% Iteration phase:
Ui = randn(nnz(I), 1); % Initial guess.

fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
if any(strcmpi(iter_type, {'MPE', 'RRE'}))
    cycle = 20; % Actually each iteration computes (cycle + 1) vectors.
    iters = iters / cycle;
    [Uf, residuals] = extrapolate(Ui, @(u) T*u + d, cycle, iters, iter_type);
else
    [Uf, residuals] = iterate(Ui, A, f, R, iters, iter_type); 
end
fprintf('(%.3fs)\n', toc);


%% Save and show the results.
figure(1); clf;
mat_file = 'results.mat';
save(mat_file)
err = show(mat_file); 
fprintf('error: %e\n', err);
