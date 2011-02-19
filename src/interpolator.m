function A = interpolator(X, Xi)
    dir = size(X) - size(Xi);
    assert(all(sort(dir) == [0 1]))
    K1 = false(size(X));
    K1(1:end-dir(1), 1:end-dir(2)) = true;
    K2 = shift(K1, dir);
    
    D = X(K2) - X(K1);
    assert(all(D(:) > 0));
    A = [(X(K2) - Xi(:)), (Xi(:) - X(K1))] ./ [D D];
    assert(all(A(:) <= 1));
    assert(all(A(:) >= 0));
    I = col(1:numel(Xi));
    A = sparse([I I], [find(K1) find(K2)], A);
end
