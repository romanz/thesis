function M = inv_diag(A, I)
    a = diag(A);
    M = sparse(I, I, 1./a(I), size(A, 1), size(A, 2));
end

