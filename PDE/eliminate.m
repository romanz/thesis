function [A, f] = eliminate(A, f, J)
% Eliminate J variables from Av = f (for boundary elimination) 
% using Gaussian Elimination.
M = [A f];
N = size(A, 1);

R1 = speye(N);
S = A(J, J);
assert(nnz(S - diag(diag(S))) == 0)
R1(J, J) = dinv(S);
M = R1 * M;

R2 = speye(N);
R2(:, J) = -A(:, J);
R2(J, J) = 0;
M = R2 * M;

A = M(:, 1:end-1);
f = M(:, end);
