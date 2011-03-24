% P = SELECT(I)
%   Select elements given by logical I, implemented by sparse matrix.
%   Can be thought as an operator of {x -> x(I)}.
function P = select(I)
    N = numel(I);
    interior = find(I);
    M = numel(interior);
    P = sparse(1:M, interior, ones(M, 1), M, N);
end
