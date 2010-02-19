%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for the solver.
x = logspace(0, 1, 1+2^5); 
y = linspace(1, 10, 1+2^5); 

% We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid(x, y);
sz = [numel(x) numel(y)];

% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It is useful for solver's verification.
U = @(X, Y) zeros(sz)  +   X + 2*Y;
L = @(X, Y) zeros(sz)  +   10;
C = @(X, Y) zeros(sz)  +   4*X + 3*Y;

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, I] = laplacian(X, Y, C(X, Y), sz); fprintf('(%.3fs)\n', toc);
V0 = U(X, Y); % The ideal solution (for boundary condition)
Vi = V0; % Initial guess.
Vi(I) = 0; % "Fill" the interior with initial guess
F = L(X, Y); % right hand side of the equation
fprintf('Apply Jacobi solver... '); tic;
Vf = jacobi(A, Vi, F, 10e3); fprintf('(%.3fs)\n', toc);

% Show the results
show(X, Y, Vf, V0)
