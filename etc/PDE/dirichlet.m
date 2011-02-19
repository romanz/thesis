function [A, f] = dirichlet(A, f, J, U)
% Substitute Dirichlet boundary conditions v(j) = u,
% into Av = f linear system, where U = @(X, Y) is a 
% function handle and J is grid variables' indices.
for k = 1:numel(J)
    j = J(k); % variable #j (to be eliminated)
    f(j) = U(k); % its value (as a boundary point)
    A(j, j) = 1; % 1 * v(j) = f(j)
end
