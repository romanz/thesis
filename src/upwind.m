function P = upwind(grid, V, dim)
    dir = (1:2 == dim);
    assert(any(dir));
    h = ones(dir + 1);
    I = grid.I;
    K1 = find(I | shift(I, -dir));
    K2 = find(I | shift(I, +dir));
    K = logical( convn(grid.I, h, 'valid') );
    % The following code maybe optimized using SPINIT scheme.
    V = convn(V, h, 'valid') < 0;
    V = V(K);
    P = sparse([1:nnz(K) 1:nnz(K)], [K1 K2], [V ~V], nnz(K), grid.numel);
end
