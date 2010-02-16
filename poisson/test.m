function test
%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used for solution.

clc
sz = 1+[2^5, 0]; % We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = create_grid( logspace(0, 1, sz(1)), linspace(-1, 1, sz(2)) );
[U, L] = create_sol(); % The ideal solution and its Laplacian.

% We actually solve the linear system: Av = f
[A, I] = create_laplacian(X, Y, sz);
V0 = U(X, Y); % The ideal solution (for boundary condition)
Vi = V0; % Initial guess.
Vi(I) = 0; % "Fill" the interior with initial guess
F = L(X, Y); % right hand side of the equation
Vf = jacobi(A, Vi, F, 10000);
% subplot 131; mesh(X, Y, V0); xlabel('X'); ylabel('Y'); 
% subplot 132; mesh(X, Y, Vi); xlabel('X'); ylabel('Y'); 
% subplot 133; mesh(X, Y, Vf-V0); xlabel('X'); ylabel('Y'); 
subplot 121; 
plot(X, [V0 Vf Vi]); 
xlabel('X'); title('Solution'); legend('V_0', 'V_f', 'V_i');
subplot 122; 
semilogy(X, abs(Vf-V0)); 
xlabel('X'); title('Error')
save poisson_test

function [v] = jacobi(A, v, f, T)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% v' = v + D^{-1} (f - Av)
d = diag(A);
I = find(d);
f = f(I); % Restrict RHS.
A = A(I, :); % Restric operator.
d = 1 ./ d(I); % Restrict the inverse.
for t = 1:T
    Av = A*v(:); % Apply A on v.
    r = f - Av; % Compute residual.
    dv = r .* d; % Compute an update.
    v(I) = v(I) + dv; % Update v's interior.
end

function [U, L] = create_sol()
%% Create solutions for the specific diff. eq. instance.
% - F is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It may be useful for solver's verification.
U = @(X, Y) (10-X).^2;
L = @(X, Y) 0*X + 2;

function [X, Y] = create_grid(x, y)
%% Create a grid specified by {X(k), Y(k)}.
% X, Y, I are 2D arrays, created by NDGRID, 
% so X is the 1st coordinate and Y is the 2nd!
[X, Y] = ndgrid(x, y);

function [A, I] = create_laplacian(X, Y, sz)
%% Laplacian discretization using sparse matrix
N = prod(sz);
I = true(sz);
if sz(1) > 1
    I(1, :) = false; 
    I(end, :) = false; 
end
if sz(2) > 1
    I(:, 1) = false; 
    I(:, end) = false; 
end

K = find(I); % Fill only interior points
A = spalloc(N, N, 5*N); % Preallocate sparse matrix.
P = [0 -1 1; 1 0 -1; -1 1 0];
for k = K(:).'
    if sz(1) > 1
        [i, j] = ind2sub(sz, k);
        i = i + [-1; 0; 1];
        j = j + [ 0; 0; 0];    
        m = sub2ind(sz, i, j);
        x = X(m);
        x = P * x; % % [x3-x2; x1-x3; x2-x1]
        x = x / prod(x);
        A(k, m) = A(k, m) + x.';
    end

    if sz(2) > 1
        [i, j] = ind2sub(sz, k);
        i = i + [ 0; 0; 0];
        j = j + [-1; 0; 1];
        m = sub2ind(sz, i, j);
        y = Y(m);
        y = P * y; % [y3-y2; y1-y3; y2-y1]
        y = y / prod(y);
        A(k, m) = A(k, m) + y.';
    end
end
A = -2*A; % Fix missing factor and sign-reverse.
