function main
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for the solver.
m = 5;
x = logspace( -1, 0, 1+2^m*2); 
y = linspace(0.1, 1, 1+2^m); 
% # of iterations
T = 12e3;
% type = 'RedBlack';
type = 'Jacobi';

% We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid(x, y);
sz = [numel(x) numel(y)];

% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It is useful for solver's verification.
U = @(X, Y) zeros(sz)  +   sin(X) + cos(Y);
C = @(X, Y) zeros(sz)  +   exp(X - Y);
L = @(X, Y) zeros(sz)  +   -exp(X - Y) .* (sin(X) - cos(X) + cos(Y) - sin(Y));

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, I] = laplacian(sz, X, Y, C(X, Y)); fprintf('(%.3fs)\n', toc);
V0 = U(X, Y); % The ideal solution (for boundary condition)
Vi = V0; % Initial guess.
Vi(I) = 0; % "Fill" the interior with initial guess

F = L(X, Y); % right hand side of the equation
fprintf('Construct iterative %s solver... ', type); tic;
[C, d] = jacobi(A, Vi, F, I); fprintf('(%.3fs)\n', toc);

fprintf('Apply %d iterations... ', T); tic;
[Vf] = iterate(Vi, C, d, T, type); fprintf('(%.3fs)\n', toc);


% Save and show the results.
mat_file = 'results.mat';
save(mat_file)
show(mat_file)
