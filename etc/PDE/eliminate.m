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

I = true(N, 1);
I(J) = false;
I = find(I);

S = A(I, J);
[r, c, v] = find(S);
R2 = sparse([I; I(r)], [I; J(c)], ...
    [ones(numel(I), 1); -v], N, N);
M = R2 * M;

A = M(:, 1:end-1);
f = M(:, end);
