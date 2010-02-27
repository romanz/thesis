function [C, d] = jacobi(A, v, f)
%% Compute {C,d} for Jacobi iteration.
% Av = [D + (A-D)]v = f
% Dv = f - (A-D)v = (f - Av) + Dv
% Dinv := D^{-1}
% v' = v + Dinv (f - Av) = (I - Dinv A)v + Dinv f = Cv + d
N = numel(v);
D = diag(A);
Dinv = sparse(1:N, 1:N, 1./D);
I = speye(N);

%   v' = Cv + d
%   C := I - Dinv * A
%   d := Dinv * f
C = I - Dinv * A;
d = Dinv * f(:);

