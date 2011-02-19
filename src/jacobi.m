function M = jacobi(sz, A)    
    M{1} = inv_diag(A, 1:prod(sz));
end
