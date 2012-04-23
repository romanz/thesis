% D = spdiag(d)
%   Create sparse diagonal matrix, using d values.
function D = spdiag(d)
n = numel(d);
d = sparse(d(:));
D = sparse(1:n, 1:n, d, n, n);