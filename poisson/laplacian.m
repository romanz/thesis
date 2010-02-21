function [A, interior] = laplacian(X, Y, C, sz)
%% Laplacian discretization using sparse matrix
N = prod(sz);
interior = true(sz);
% Compute interior points' indices (special handling of 1D)
if sz(1) > 1
    interior(1, :) = false; 
    interior(end, :) = false; 
end
if sz(2) > 1
    interior(:, 1) = false; 
    interior(:, end) = false; 
end
ind = @(I, J) sub2ind(sz, I, J);
K = find(interior); % Fill only interior points
K = K(:);
[I, J] = ind2sub(sz, K);

Kr = ind(I+1, J);
Kl = ind(I-1, J);
Ku = ind(I, J+1);
Kd = ind(I, J-1);

% Build Laplacian sparse matrix:
Dxx = ...
 ( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ).*( X(Kr) - X(Kl) )) * [1 -1 0] -  ...
 ( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ).*( X(Kr) - X(Kl) )) * [0 1 -1];

Dyy = ...
 ( C(Ku) + C(K) ) ./ (( Y(Ku) - Y(K) ).*( Y(Ku) - Y(Kd) )) * [1 -1 0] -  ...
 ( C(K) + C(Kd) ) ./ (( Y(K) - Y(Kd) ).*( Y(Ku) - Y(Kd) )) * [0 1 -1];


D = repmat([1 0 -1], [numel(K) 1]);
I = repmat(I, [1 3]);
J = repmat(J, [1 3]);
K = repmat(K, [1 3]);

P = ind(I + D, J);
Dxx = sparse(K, P, Dxx, N, N);

P = ind(I, J + D);
Dyy = sparse(K, P, Dyy, N, N);

A = Dxx + Dyy;
% full(A) 
