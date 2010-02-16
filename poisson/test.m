%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used for solution.
clc

% Create solutions for the specific diff. eq. instance.
% - F is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It may be useful for solver's verification.
U = @(X, Y) (10-X).^2;
L = @(X, Y) 0*X + 2;

% Create grid for solver
sz = 1+[2^5, 0]; % We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid( logspace(0, 1, sz(1)), linspace(-1, 1, sz(2)) );

% We actually solve the linear system: Av = f
[A, I] = laplacian(X, Y, sz);
V0 = U(X, Y); % The ideal solution (for boundary condition)
Vi = V0; % Initial guess.
Vi(I) = 0; % "Fill" the interior with initial guess
F = L(X, Y); % right hand side of the equation
Vf = jacobi(A, Vi, F, 10000);

% Show results
subplot 121; 
plot(X, [V0 Vf Vi]); 
xlabel('X'); title('Solution'); legend('V_0', 'V_f', 'V_i');
subplot 122; 
semilogy(X, abs(Vf-V0)); 
xlabel('X'); title('Error')
