function V = stagger(sz, X, Y, dim, func)
assert(sz(dim) > 1)

K = 1:numel(sz) == dim;
h = ones(K + 1);
h = h / sum(h); % Averaging filter

X = average(X, h);
Y = average(Y, h);
V = func(X, Y);

switch dim
    case 1
        V = V(:, 2:end-1);
    case 2
        V = V(2:end-1, :);
    otherwise
        error('2D support only!')
end
