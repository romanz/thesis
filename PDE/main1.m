function main1
fprintf('\n');
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

%% Create grid for the solver.
m = 5;
x = linspace(0, 1, 1+2^m); 
y = linspace(0, 1, 1+2^m); 

%% # of iterations
iters = 10e3;

iter_type = '';
iter_type = 'Jacobi';
% iter_type = 'RedBlack'; iters = iters / 2;
% iter_type = 'RRE'; iters = 4e3;
% iter_type = 'MPE'; iters = 4e3;

% Dirichlet/Nemunann
conditions.left  = 'Dirichlet';
conditions.right = 'Dirichlet';
conditions.up    = 'Dirichlet'; % 'Neumann'; 
conditions.down  = 'Dirichlet';

%% We use NDGRID convention (X is 1st, Y is 2nd)
sz = [numel(x) numel(y)];
N = prod(sz);
Z = zeros(sz);
I = interior(sz);
[X, Y] = ndgrid(x, y);
Vx = stagger(I, X, Y, 1, @(x, y) x+y);
Vy = stagger(I, X, Y, 2, @(x, y) x-y);
alpha = 1;
%% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the diverence of C * grad(U).
% It is useful for solver's verification.
U = @(X, Y) X .* Y;
Ux = @(X, Y) Y;
Uy = @(X, Y) X;
C = @(X, Y) 1 + X .* Y;
F = @(X, Y) Z;

%% We actually solve the linear system: Av = f
fprintf('Compute Laplacian on %d x %d grid... ', numel(x), numel(y)); tic;
diffusion = laplacian(I, X, Y, C(X, Y));
fi = F(X, Y); % Right-hand side of the equation
Ai = diffusion - alpha * advection(I, X, Y, Vx, Vy, 'central');
fi = fi(:); 
fi(~I) = 0;

Bl = boundary(I, [-1 0]); Il = shift(Bl, [+1 0]);
Br = boundary(I, [+1 0]); Ir = shift(Br, [-1 0]);
Bd = boundary(I, [0 -1]); Id = shift(Bd, [0 +1]);
Bu = boundary(I, [0 +1]); Iu = shift(Bu, [0 -1]);

% Add boundary conditions for X:
Ab = spalloc(N, N, N);
fb = zeros(N, 1);
if sz(1) > 1
    switch conditions.left
        case 'Dirichlet'
            [Ab, fb] = dirichlet( Ab, fb, find(Bl), U(X(Bl), Y(Bl)) );
        case 'Neumann'
            [Ab, fb] = neumann(Ab, fb, find(Bl), find(Il), ...
                Ux((X(Bl) + X(Il))/2, (Y(Bl) + Y(Il))/2) .* (X(Bl) - X(Il)));
    end
    switch conditions.right
        case 'Dirichlet'
            [Ab, fb] = dirichlet( Ab, fb, find(Br), U(X(Br), Y(Br)) );
        case 'Neumann'
            [Ab, fb] = neumann(Ab, fb, find(Br), find(Ir), ...
                Ux((X(Br) + X(Ir))/2, (Y(Br) + Y(Ir))/2) .* (X(Br) - X(Ir)));
    end
end
% Add boundary conditions for Y:
if sz(2) > 1
    switch conditions.down
        case 'Dirichlet'
            [Ab, fb] = dirichlet( Ab, fb, find(Bd), U(X(Bd), Y(Bd)) );
        case 'Neumann'
            [Ab, fb] = neumann(Ab, fb, find(Bd), find(Id), ...
                Uy((X(Bd) + X(Id))/2, (Y(Bd) + Y(Id))/2) .* (Y(Bd) - Y(Id)));
    end
    switch conditions.up
        case 'Dirichlet'
            [Ab, fb] = dirichlet( Ab, fb, find(Bu), U(X(Bu), Y(Bu)) );
        case 'Neumann'
            [Ab, fb] = neumann(Ab, fb, find(Bu), find(Iu), ...
                Uy((X(Bu) + X(Iu))/2, (Y(Bu) + Y(Iu))/2) .* (Y(Bu) - Y(Iu)));
    end
end
assert(nnz(Ai & Ab) == 0);
assert(nnz(fi & fb) == 0);
% Eliminate boundary variables, forming linear system only for the interior.
[Ae, fe] = eliminate(Ai + Ab, fi + fb, find(~I));
[Ar, fr] = restrict(Ae, fe, I, false);
fprintf('(%.3fs)\n', toc);

fprintf('Construct Jacobi iteration... '); tic;
[R, T, d] = jacobi(Ar, fr); fprintf('(%.3fs)\n', toc);

% fprintf('Maximal eigenvalue of T: ');
% lambda = abs(max_eigs(T, 1));
% fprintf('%.10f => %d iterations/decade\n', ...
%     lambda, ceil(-1/log10(lambda)));

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
    Uf = Ar \ fr;
end
iter_time = toc;
fprintf('(%.3fs) [%.1fns/cell]\n', ...
    iter_time, 1e9*iter_time/(iters*N)); % ~50ns

% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
err = show(load(mat_file));
fprintf('error: %e\n', err);

%% Expected output:
% Compute Laplacian on 65 x 33 grid... (2.747s)
% Construct Jacobi iteration... (0.166s)
% Maximal eigenvalue of T: 0.9988276216 => 1963 iterations/decade
% Apply Jacobi solver [30000]... (5.132s)
% error: 8.374371e-006
