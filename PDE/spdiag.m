function D = spdiag(X)

X = sparse(X(:));
D = diag(X);
