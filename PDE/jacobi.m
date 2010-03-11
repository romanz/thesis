function [R, T, d] = jacobi(A, f) %#ok<INUSD>
%% Compute the appropriate matrices for for Jacobi iteration.
% Av = [D + (A-D)]v = f
% Dv = f - (A-D)v = (f - Av) + Dv
% R := D^{-1}
% v' = v + R (f - Av) = (I - R A)v + R f
% v' = Tv + d
d = diag(A); % Vector of diagonal elements of A
N = numel(d);
R = sparse(1:N, 1:N, 1./d);
T = speye(N) - R * A;
d = R * f;
