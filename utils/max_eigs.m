function e = max_eigs(A, k)
% Compute maximal k eigenvalues of sparse matrix A.
if nnz(A) > 0
    S = struct('disp', 0, 'v0', ones(length(A), 1));
    e = eigs(A, k, 'lm', S); 
else
    e = 0;
end
