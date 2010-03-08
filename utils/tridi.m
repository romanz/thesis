function A = tridi(L, D, U)
N = numel(D);
A = sparse(N, N);
A(sub2ind(size(A), 1:N, 1:N)) = D;
A(sub2ind(size(A), 1:N-1, 2:N)) = U;
A(sub2ind(size(A), 2:N, 1:N-1)) = L;
