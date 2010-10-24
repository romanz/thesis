% P = EXPAND(I)
%   Expand interior by zero padding, implemented by sparse matrix.
%   Can be thought as an adjoint operator of {x -> x(I)}.
function P = expand(I)
    N = numel(I);
    interior = find(I);
    M = numel(interior);
    P = sparse(interior, 1:M, ones(M, 1), N, M);
end
