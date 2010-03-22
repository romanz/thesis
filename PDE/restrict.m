function [A, f] = restrict(A, f, I, self_adjoint)
%% Sanity checks for the linear system
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0); 
assert(nnz(A(:, ~I)) == 0); 
assert(nnz(f(~I)) == 0);
% Verify that all entries are valid real numbers:
assert(nnz(isnan(A) | isinf(A)) == 0)
% Verify that the matrix is symmetric (the original operator is self-adjoint):
if nargin >= 4 && self_adjoint
    assert(nnz(A - A') == 0)
end
%% Restrict the problem to interior variables
A = A(I, I); 
f = f(I);
