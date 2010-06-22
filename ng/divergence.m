% Sparse divergence operator.
function G = divergence(K, X, Y, dir)
    sz = size(K);
    [I, J] = ind2sub(sz, find(K));
    K1 = sub2ind(sz, I, J);
    K2 = sub2ind(sz, I + dir(1), J + dir(2));
    K0 = sub2ind(sz - dir, I, J);
    Xc = (X(K1) + X(K2))/2;
    Yc = (Y(K1) + Y(K2))/2;
    switch find(dir)
        case 1
            A1 = X(K1) .^ 2;
            A2 = X(K2) .^ 2;
            D = (X(K2) - X(K1)) .* (Xc .^ 2);
        case 2
            A1 = sin(Y(K1));
            A2 = sin(Y(K2));
            D = (Y(K2) - Y(K1)) .* (Xc .* sin(Yc));
    end
    G = sparse([1:numel(K0), 1:numel(K0)], [K1 K2], [-A1./D,  A2./D], ...
        numel(K0), prod(sz));
end
