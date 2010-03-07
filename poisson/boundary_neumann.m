function [A, F] = boundary_neumann(A, F, J, delta, X, Y, Ux, Uy)
sz = size(J);
J = find(J);
for k = 1:numel(J)
    % Jb - boundary variable (to be eliminated)
    % Ji - interior variable, such that: U(Jb) - U(Ji) = dU
    Jb = J(k);
    Ji = shift_index(sz, Jb, -delta);
    Xc = (X(Jb) + X(Ji))/2;
    Yc = (Y(Jb) + Y(Ji))/2;
    dU = Ux(Xc, Yc) .* (X(Jb) - X(Ji)) + ...
         Uy(Xc, Yc) .* (Y(Jb) - Y(Ji));
    I = find( A(:, Jb) ); % corresponding equations
    F(I) = F(I) - A(I, Jb) * dU; % update RHS
    A(I, Ji) = A(I, Ji) + A(I, Jb); % update LHS
    A(I, Jb) = 0; % Eliminate this variable from A.
end

function index = shift_index(sz, index, delta)
[i, j] = ind2sub(sz, index);
di = delta(1);
dj = delta(2);
index = sub2ind(sz, i + di, j + dj);
