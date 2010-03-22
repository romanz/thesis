function [A, f] = dirichlet(A, f, J, U)
% Substitute Dirichlet boundary conditions v(j) = u,
% into Av = f linear system, where U = @(X, Y) is a 
% function handle and J is grid variables' indices.
for k = 1:numel(J)
    j = J(k); % variable #j (to be eliminated)
    u = U(k); % its value (as a boundary point)
    I = find(A(:, j)); % corresponding equations
    f(I) = f(I) - A(I, j) * u; % eliminate and update RHS
    A(I, j) = 0; % Update A (this variable is unused)
end
