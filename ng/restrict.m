% Keep only the columns indexed by I (logical mask)
function P = restrict(I)
    N = numel(I);
    J = find(I);
    P = sparse(J, 1:numel(J), ones(size(J)), N, numel(J));
end