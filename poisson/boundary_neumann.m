function [A, F] = boundary_neumann(A, F, Jb, Ji, dV)
if islogical(Jb)
    Jb = find(Jb);
end
if islogical(Ji)
    Ji = find(Ji);
end
for k = 1:numel(dV)
    % Jb(k) - boundary variable (to be eliminated)
    % Ji(k) - interior variable, such that: V(Jb(k)) - V(Ji(k)) = dV(k)
    I = find( A(:, Jb(k)) ); % corresponding equations
    F(I) = F(I) - A(I, Jb(k)) * dV(k); % update RHS
    A(I, Ji(k)) = A(I, Ji(k)) + A(I, Jb(k)); % update LHS
    A(I, Jb(k)) = 0; % Eliminate this variable from A.
end
