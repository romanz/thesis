function [R, T, d] = jacobi(A, f) %#ok<INUSD>
%% Compute the appropriate matrices for for Jacobi iteration.
% Av = [D + (A-D)]v = f
% Dv = f - (A-D)v = (f - Av) + Dv
% R := D^{-1}
% v' = v + R (f - Av) = (I - R A)v + R f
% v' = Tv + d
N = numel(f);
D = diag(A);
R = sparse(1:N, 1:N, 1./D);
T = speye(N) - R * A;
d = R * f;
