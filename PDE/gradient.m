% Sparse gradient operator.
function G = gradient(sz, X, Y, dim)

K = true(sz);
switch dim
    case 1, K(end, :) = false; 
    case 2, K(:, end) = false; 
end
[I, J] = ind2sub(sz, find(K));
P = (1:numel(sz)) == dim;
K1 = sub2ind(sz, I, J);
K2 = sub2ind(sz, I + P(1), J + P(2));
K0 = sub2ind(sz - P, I, J);
switch dim
    case 1, D = X(K2) - X(K1);
    case 2, D = Y(K2) - Y(K1);
end
G = sparse([K0 K0], [K1 K2], 1./[-D D]);