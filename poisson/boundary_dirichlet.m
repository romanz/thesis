function [A, F] = boundary_dirichlet(A, F, J, X, Y, U)
if islogical(J)
    J = find(J);
end
for k = 1:numel(J)
    j = J(k); % variable #j (to be eliminated)
    u = U(X(j), Y(j)); % its value (as a boundary point)
    I = find(A(:, j)); % corresponding equations
    F(I) = F(I) - A(I, j) * u; % eliminate and update RHS
    A(I, j) = 0; % Update A (this variable is unused)
end
