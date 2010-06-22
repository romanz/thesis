% Sparse gradient operator.
function G = gradient(K, X, Y, dir)
    sz = size(K);
    [I, J] = ind2sub(sz, find(K));
    K1 = sub2ind(sz, I, J);
    K2 = sub2ind(sz, I + dir(1), J + dir(2));
    K0 = sub2ind(sz - dir, I, J);
    switch find(dir)
        case 1, D = X(K2) - X(K1);
        case 2, D = (Y(K2) - Y(K1)) .* (X(K1) + X(K2))/2;
    end
    G = sparse([1:numel(K0), 1:numel(K0)], [K1 K2], [-1./D,  1./D], ...
        numel(K0), prod(sz));
end
