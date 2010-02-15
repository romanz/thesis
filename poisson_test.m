function poisson_test
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: ?u = f
% ? is discretized on a grid, and Jacobi iteration is used for solution.
clc
sz = [3, 1]; % We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = create_grid( linspace(0, 1, sz(1)), linspace(0, 1, sz(2)) );
[U, L] = create_sol(); % The ideal solution and its Laplacian.

% A is discrete version of ?, so we actually solve the system: Av = f
[A, I] = create_laplacian(X, Y, sz);
V0 = U(X, Y); % The ideal solution (for boundary condition)
V0(I) = 0; % "Fill" the interior with initial guess
F = L(X, Y); % right hand side of the equation
V = jacobi(A, V0, F, I, 1);

save poisson_test

function [v] = jacobi(A, v, f, I, T)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% v' = v + D^{-1} (f - Av)
d = diag(A);
full(d)
v
for t = 1:T
    Av = A*v; % Apply ? on v.
    r = f - Av; % Compute residual.
    dv = r(I) ./ d(I); % Compute an update.
    v(I) = v(I) + dv; % Update v's interior.
end
v

function [A, I] = create_laplacian(X, Y, sz)
%% Laplacian discretization using sparse matrix
N = prod(sz);
A = sparse(N, N);
K = true(sz);
if sz(1) > 1, K(1, :) = false; K(end, :) = false; end
if sz(2) > 1, K(:, 1) = false; K(:, end) = false; end
K = find(K); % Fill only interior points' corresponding rows
for k = K(:).'
    [kx, ky] = ind2sub(sz, k);
%     kx = kx + [-1 0 0 0 1];
%     ky = ky + [0 -1 0 1 0];
%     A(k, sub2ind(sz, kx, ky)) = ...
%               [1 1 -4 1 1];
    kx = kx + [-1 0 1];
    ky = ky + [0 0 0];    
    A(k, sub2ind(sz, kx, ky)) = ...
              [1 -2 1];
end
I = any(A, 2); % Find out interior points

function [X, Y] = create_grid(x, y)
%% PDE utility
% Create a grid specified by {X(k), Y(k)}.
% X, Y, I are 2D arrays, created by NDGRID, 
% so X is the 1st coordinate and Y is the 2nd!
[X, Y] = ndgrid(x, y);

function [U, L] = create_sol()
%% Create solutions for the specific diff. eq. instance.
% - F is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It may be useful for solver's verification.
a = 2;
U = @(X, Y) a*X+1;
L = @(X, Y) 0*X;

