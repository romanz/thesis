function [v] = jacobi(A, v, f, T)
%% Jacobi iteration:
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% v' = v + D^{-1} (f - Av)
d = diag(A);
I = find(d);
f = f(I); % Restrict RHS.
A = A(I, :); % Restric operator.
d = 1 ./ d(I); % Restrict the inverse.

% Keep the original size and convert to column vectors:
sz = size(v); 
v = v(:);
f = f(:);
for t = 1:T
    Av = A*v; % Apply A on v.
    r = f - Av; % Compute residual.
    dv = r .* d; % Compute an update.
    v(I) = v(I) + dv; % Update v's interior.
end
v = reshape(v, sz);
