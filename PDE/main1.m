function main1
fprintf('\n');
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

%% Create grid for the solver.
m = 5;
x = logspace(-1, 0, 1+2*2^m); 
y = 1+linspace(-1, 0, 1+2^m); 

%% # of iterations
iters = 30e3;

iter_type = '';
iter_type = 'Jacobi';
% iter_type = 'RedBlack'; iters = iters / 2;
% iter_type = 'RRE'; iters = 4e3;
% iter_type = 'MPE'; iters = 4e3;
conditions.left  = 'Dirichlet';
conditions.right = 'Dirichlet';
conditions.up    = 'Neumann';
conditions.down  = 'Neumann';


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
[A, M, I] = laplacian(sz, X, Y, C(X, Y));
f = L(X, Y); % Right-hand side of the equation
f = M * f(:); % Multiply by preconditioner

Bl = boundary(I, [-1 0]); Il = shift(Bl, [+1 0]);
Br = boundary(I, [+1 0]); Ir = shift(Br, [-1 0]);
Bd = boundary(I, [0 -1]); Id = shift(Bd, [0 +1]);
Bu = boundary(I, [0 +1]); Iu = shift(Bu, [0 -1]);

% Add boundary conditions for X:
if sz(1) > 1
    switch conditions.left
        case 'Dirichlet'
            [A, f] = dirichlet( A, f, find(Bl), U(X(Bl), Y(Bl)) );
        case 'Neumann'
            [A, f] = neumann(A, f, find(Bl), find(Il), ...
                Ux((X(Bl) + X(Il))/2, (Y(Bl) + Y(Il))/2) .* (X(Bl) - X(Il)));
    end
    switch conditions.right
        case 'Dirichlet'
            [A, f] = dirichlet( A, f, find(Br), U(X(Br), Y(Br)) );
        case 'Neumann'
            [A, f] = neumann(A, f, find(Br), find(Ir), ...
                Ux((X(Br) + X(Ir))/2, (Y(Br) + Y(Ir))/2) .* (X(Br) - X(Ir)));
    end
end
% Add boundary conditions for Y:
if sz(2) > 1
    switch conditions.down
        case 'Dirichlet'
            [A, f] = dirichlet( A, f, find(Bd), U(X(Bd), Y(Bd)) );
        case 'Neumann'
            [A, f] = neumann(A, f, find(Bd), find(Id), ...
                Uy((X(Bd) + X(Id))/2, (Y(Bd) + Y(Id))/2) .* (Y(Bd) - Y(Id)));
    end

    switch conditions.up
        case 'Dirichlet'
            [A, f] = dirichlet( A, f, find(Bu), U(X(Bu), Y(Bu)) );
        case 'Neumann'
            [A, f] = neumann(A, f, find(Bu), find(Iu), ...
                Uy((X(Bu) + X(Iu))/2, (Y(Bu) + Y(Iu))/2) .* (Y(Bu) - Y(Iu)));
    end

    
end
% Eliminate boundary variables, forming linear system only for the interior.
[A, f] = eliminate(A, f, find(~I));

%% Sanity checks for the linear system
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(nnz(f(~I)) == 0);
% Verify that all entries are valid real numbers:
assert(nnz(isnan(A) | isinf(A)) == 0)
% Verify that the matrix is symmetric (the original operator is self-adjoint):
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
Ui = randn(nnz(I), 1); % Initial guess.

fprintf('Apply %s solver [%d]... ', iter_type, iters); tic;
if any(strcmpi(iter_type, {'MPE', 'RRE'}))
    cycle = 20; % Actually each iteration computes (cycle + 1) vectors.
    iters = iters / cycle;
    [Uf, residuals] = extrapolate(Ui, @(u) T*u + d, cycle, iters, iter_type);
elseif ~isempty(iter_type)
    [Uf, residuals] = iterate(Ui, T, d, iters, iter_type); 
else
    Uf = A \ f;    
end
fprintf('(%.3fs)\n', toc);

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
err = show(load(mat_file));
fprintf('error: %e\n', err);
% full(A), f
