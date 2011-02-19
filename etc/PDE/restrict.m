function [A, f] = restrict(A, f, I, self_adjoint)
%% Sanity checks for the linear system
% Verify that boundary variables are eliminated properly.
assert(nnz(A(~I, :)) == 0, 'Boundary equations are not eliminated!'); 
assert(nnz(A(:, ~I)) == 0, 'Boundary variables are not eliminated!'); 
assert(nnz(f(~I)) == 0, 'Bundary RHS is not eliminated!');
% Verify that all entries are valid real numbers:
assert(nnz(isnan(A) | isinf(A)) == 0, 'NaN/inf is found in A!')
% Verify that the matrix is symmetric (the original operator is self-adjoint):
if nargin >= 4 && self_adjoint
    assert(nnz(A - A') == 0, 'A is not self-adjoint!')
end
%% Restrict the problem to interior variables
A = A(I, I); 
f = f(I);
