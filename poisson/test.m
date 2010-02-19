%% Simulate and solve Poisson Equation on a 2D grid using Relaxation
% The equation to solve is: Laplacian(u) = f.
% Laplacian is discretized on a grid, and Jacobi iteration is used.

% Create grid for solver
sz = 1+[2.^4, 2.^4]; % We use NDGRID convention (X is 1st, Y is 2nd)
[X, Y] = ndgrid( linspace(-1, 1, sz(1)), linspace(-1, 1, sz(2)) );

% Create solutions for the specific diff. eq. instance.
% - U is the function itself (for boundary conditions).
% - L is the Laplacian of F.
% It is useful for solver's verification.
U = @(X, Y) 0*X  +   X.^3 .* Y - X .* Y.^3 + 0.5;
L = @(X, Y) 0*X;

% We actually solve the linear system: Av = f
fprintf('Compute Laplacian operator... '); tic;
[A, I] = laplacian(X, Y, sz); fprintf('(%.3fs)\n', toc);
V0 = U(X, Y); % The ideal solution (for boundary condition)
Vi = V0; % Initial guess.
Vi(I) = 0; % "Fill" the interior with initial guess
F = L(X, Y); % right hand side of the equation
fprintf('Apply Jacobi solver... '); tic;
Vf = jacobi(A, Vi, F, 5000); fprintf('(%.3fs)\n', toc);

%% Show results
subplot 121; 
if sz(1) > 1 && sz(2) > 1
    mesh(X, Y, Vf); xlabel('X'); ylabel('Y'); title('Solution'); 
else
    if sz(1) > 1
        plot(X, Vf(:));
    else
        plot(Y, Vf(:));
    end
end

subplot 122; 
E = Vf - V0; norm(E, inf)
if sz(1) > 1 && sz(2) > 1
    mesh(X, Y, abs(E));  xlabel('X'); ylabel('Y'); title('Error')
else
    if sz(1) > 1
        plot(X, abs(E(:)));
    else
        plot(Y, abs(E(:)));
    end
end

