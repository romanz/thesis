% Expand interior by zero padding.
function P = expand(I)
    N = numel(I);
    interior = find(I);
    M = numel(interior);
    P = sparse(interior, 1:M, repmat(1, [M 1]), N, M);
end
