function [P, V] = upwind(grid, V, dim)
    dir = (1:2 == dim);
    assert(any(dir));
    h = ones(dir + 1);
    q = ~dir; % remove other dimension's velocity ghost points
    I = grid.I(1+q(1):end-q(1), 1+q(2):end-q(2));
    K1 = find(I | shift(I, -dir));
    K2 = find(I | shift(I, +dir));
    K = logical( convn(I, h, 'valid') );
    % The following code maybe optimized using SPINIT scheme.
    Y = convn(V, h, 'valid') < 0;
    Y = Y(K);
    P = sparse([1:nnz(K) 1:nnz(K)], [K1 K2], [Y ~Y], nnz(K), numel(I));
end
