function M = redblack(sz, A)    
    K = checkerboard(sz);
    I0 = find(K == 0);
    I1 = find(K == 1);
    M{1} = inv_diag(A, I0);
    M{2} = inv_diag(A, I1);
end

function M = inv_diag(A, I)
    a = diag(A);
    M = sparse(I, I, 1./a(I), size(A, 1), size(A, 2));
end
