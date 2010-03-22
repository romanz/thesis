fprintf('\n'); 
clear;
%% Simulate and solve Diffusion-Advection Equation on a 2D grid.
% The equation to solve is: D Laplacian(u) - V*Gradient(u) = f.
% Central difference is used for Laplacian. 
% Upstream difference is used for Gradient.
% Jacobi iteration is used for relaxation.
randn('state', 1);
iters = 10e3;
iter_type = 'RedBlack';
sz = 1+2.^([1 1]*5);

%% Create the grid
N = prod(sz);
x = 2*linspace(-1, 1, sz(1));
y = 2*linspace(-1, 1, sz(2));
[X, Y] = ndgrid(x, y);

%% The discrete Laplacian operator
[L, M, I] = laplacian(sz, X, Y);
Diffusion = dinv(M) * L;

%% Velocity field (Vx, Vy) definition on a staggered grid
Vx = average( Y, [1; 1]/2); % Evaluated at left and right edges.
Vy = average(-X, [1, 1]/2); % Evaluated at top and bottom edges.
if sz(2) > 1 && numel(Vx)
    Vx = average(Vx, [0, 1, 0]); % Remove Y boundary
end
if sz(1) > 1 && numel(Vy)
    Vy = average(Vy, [0; 1; 0]); % Remove X boundary
end

%% The discrete advection operator
Advection = gradient(sz, X, Y, Vx, Vy, 'centered');

%% The operator combined operator
A = Diffusion - Advection; 
U = @(X, Y) exp(-(X.^2 + Y.^2) / 2);
f = (X.^2 + Y.^2 - 2) .* exp(-(X.^2 + Y.^2) / 2); % Right-Hand Side
f(~I) = NaN; % Not defined on the boundary

%% Substitute Boundary Conditions
[A, f] = dirichlet(A, f, I, [-1 0], X, Y, U);
[A, f] = dirichlet(A, f, I, [+1 0], X, Y, U);
[A, f] = dirichlet(A, f, I, [0 -1], X, Y, U);
[A, f] = dirichlet(A, f, I, [0 +1], X, Y, U);
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(all(isnan(f(~I))) == 1);
assert(nnz(isnan(A(I, I))) == 0);
assert(nnz(isnan(f(I))) == 0);
% Restrict the operator and RHS on the interior variables
A = A(I, I); 
f = f(I);

%% Solve the linear system
if isempty(iter_type) % exact solver
    Uf = A \ f(:);
else
% Apply relaxation to solve the problem
    fprintf('Construct Jacobi iteration... '); tic;
    [R, T, d] = jacobi(A, f); fprintf('(%.3fs)\n', toc);

    % Iteration phase:
    Ui = randn(nnz(I), 1); % Initial guess.

    fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
    if any(strcmpi(iter_type, {'MPE', 'RRE'}))
        cycle = 20; % Actually each iteration computes (cycle + 1) vectors.
        iters = iters / cycle;
        [Uf, residuals] = extrapolate(Ui, @(u) T*u + d, cycle, iters, iter_type);
    else
        [Uf, residuals] = iterate(Ui, T, d, iters, iter_type); 
    end
end
fprintf('(%.3fs)\n', toc);

%% Save and show the results.
figure(1); clf;
mat_file = 'results.mat';
delete(mat_file);
save(mat_file);
err = show(load(mat_file)); 
fprintf('error: %e\n', err);
