function [A, F] = boundary_dirichlet(A, F, J, X, Y, U)
% Substitute Dirichlet boundary conditions v(j) = U(X(j), Y(j)) 
% into Av = f linear system, where U = @(X, Y) is a function handle
% and J is logical boundary grid.
J = find(J);
for k = 1:numel(J)
    j = J(k); % variable #j (to be eliminated)
    u = U(X(j), Y(j)); % its value (as a boundary point)
    I = find(A(:, j)); % corresponding equations
    F(I) = F(I) - A(I, j) * u; % eliminate and update RHS
    A(I, j) = 0; % Update A (this variable is unused)
end
