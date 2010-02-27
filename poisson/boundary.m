function [A, F] = boundary(A, F, J, V)
if islogical(J)
    J = find(J);
end
for k = 1:numel(J)
    j = J(k); % variable #j (to be eliminated)
    v = V(k); % its value (as a boundary point)
    I = find(A(:, j)); % corresponding equations
    F(I) = F(I) - A(I, j) * v; % eliminate and update RHS
    A(I, j) = 0; % Update A (this variable is unused)
end
