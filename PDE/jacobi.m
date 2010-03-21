function [R, T, d] = jacobi(A, f) %#ok<INUSD>
%% Compute the appropriate matrices for for Jacobi iteration.
% Au = [D + (A-D)]u = f
% Du = f - (A-D)u = (f - Au) + Du
% R := D^{-1}
% u' = u + R (f - Au) = (I - R A)u + R f
% u' = Tu + d
d = diag(A); % Vector of diagonal elements of A
N = numel(d);
R = sparse(1:N, 1:N, 1./d);
T = speye(N) - R * A;
d = R * f(:);
