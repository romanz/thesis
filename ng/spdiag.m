function D = spdiag(d)
n = numel(d);
D = sparse(1:n, 1:n, d, n, n);