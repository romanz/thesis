function [v] = jacobi(A, v, f, T)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% v' = v + D^{-1} (f - Av)
I = any(A, 2); % the interior of the grid (has non-zero on its row).
d = diag(A);

f = f(I); % Restrict RHS.
A = A(I, :); % Restric operator.

k = nnz(I); % # of the interior points in the grid.
R = sparse(1:k, 1:k, 1 ./ d(I)); % Restrict the inverse.

% Keep the original size and convert to column vectors:
sz = size(v); 
v = v(:);
f = f(:);
for t = 1:T
    Av = A*v; % Apply A on v.
    r = f - Av; % Compute residual.
    dv = R * r; % Compute an update.
    v(I) = v(I) + dv; % Update v's interior.
end
v = reshape(v, sz);
