function [v] = jacobi(A, v, f, interior, T)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% J := D^{-1}
% v' = v + J (f - Av) = (I - JA)v + Jf = Rv + d
D = diag(A);
f = f(interior); % Restrict RHS.
A = A(interior, :); % Restric operator.

k = nnz(interior);
J = sparse(1:k, 1:k, 1 ./ D(interior)); % Restrict the inverse.

R = speye(size(A, 2));
R = R(interior, :) - J * A;
d = J * f;

% Keep the original size and convert to column vectors:
sz = size(v); 
v = v(:);
f = f(:);
for t = 1:T
    v(interior) = R * v + d; % Update v's interior.
end
v = reshape(v, sz);
