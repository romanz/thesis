function R = dinv(A)
% R = pseudo-inverse of diag(A).
d = diag(A);
I = find(d);
d(I) = 1 ./ d(I);
N = numel(d);
R = sparse(1:N, 1:N, d);
