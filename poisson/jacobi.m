function [C, d] = jacobi(A, v, f, interior)
%% Compute {C,d} for Jacobi iteration.
% Av = (D+T)v = f
% Dv = f-Tv = (f - Av) + Dv
% J := D^{-1}
% v' = v + J (f - Av) = (I - JA)v + Jf = Cv + d

N = numel(v);
D = diag(A);
f = f(interior); % Restrict RHS.
A = A(interior, :); % Restric operator.

k = nnz(interior);
J = sparse(1:k, 1:k, 1 ./ D(interior)); % Restrict the inverse.
I = speye(N);

% Interior update:
%   v' = Cv + d
%   C := I - JA (on the interior)
%   d := Jf
C_i = speye(N); % Keep boundary values unchanged.
d_i = zeros(N, 1);
C_i(interior, :) = I(interior, :) - J * A;
d_i(interior) = J * f(:);

% Boundary update:
%   v' = Cv + d
%   C := 0 (on the boundary)
%   d := (Dirichlet boundary conditions)
C_b = speye(N); % Keep interior values unchanged.
C_b(~interior, :) = 0;
d_b = zeros(N, 1);
d_b(~interior) = v(~interior);

% Composition of a boundary update, followed by an interior update
% v' = Ci(Cb v + db) + di = C_i C_b v + C_i d_b + d_i
C = C_i * C_b;
d =  C_i * d_b + d_i;

