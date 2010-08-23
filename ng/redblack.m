function M = redblack(sz, A)    
    K = checkerboard(sz);
    I0 = find(K == 0);
    I1 = find(K == 1);
    M{1} = inv_diag(A, I0);
    M{2} = inv_diag(A, I1);
end
