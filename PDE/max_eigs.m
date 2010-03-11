function e = max_eigs(A, k)
% Compute maximal k eigenvalues of sparse matrix A.
if nnz(A) > 0
    e = eigs(A, k, 'lm', struct('disp', 0)); 
else
    e = 0;
end
