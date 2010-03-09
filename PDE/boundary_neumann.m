function [A, F] = boundary_neumann(A, F, J, P, X, Y, Ux, Uy)
% Substitute Neumann boundary conditions into Av = f linear system, 
% where (Ux, Uy) is the gradient, J is one of the grid's boundaries and 
% P is an unit vector, normal to the boundary, pointing outside of the grid.
sz = size(J);
J = find(J);
for k = 1:numel(J)
    % Jb - boundary variable (to be eliminated)
    % Ji - interior variable, closest to Jb.
    Jb = J(k);
    Ji = shift_index(sz, Jb, -P);
    % (Xc, Yc) is the midpoint between boundary and the interior points.
    Xc = (X(Jb) + X(Ji))/2;
    Yc = (Y(Jb) + Y(Ji))/2;
    % We use gradient at (Xc, Yc) to approximate dU = U(Jb) - U(Ji):
    dU = Ux(Xc, Yc) .* (X(Jb) - X(Ji)) + ...
         Uy(Xc, Yc) .* (Y(Jb) - Y(Ji));
    % NOTE: this approximation is accurate for quadratic U.    
    I = find( A(:, Jb) ); % corresponding equations for U(Jb).
    F(I) = F(I) - A(I, Jb) * dU; % update RHS
    A(I, Ji) = A(I, Ji) + A(I, Jb); % update LHS
    A(I, Jb) = 0; % Eliminate this variable from A.
end

function index = shift_index(sz, index, delta)
[i, j] = ind2sub(sz, index);
di = delta(1);
dj = delta(2);
index = sub2ind(sz, i + di, j + dj);
