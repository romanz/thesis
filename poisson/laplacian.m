function [A, I] = laplacian(X, Y, C, sz)
%% Laplacian discretization using sparse matrix
N = prod(sz);
I = true(sz);
% Compute interior points' indices (special handling of 1D)
if sz(1) > 1
    I(1, :) = false; 
    I(end, :) = false; 
end
if sz(2) > 1
    I(:, 1) = false; 
    I(:, end) = false; 
end

K = find(I); % Fill only interior points
% Preallocate sparse matrices:
Dxx = spalloc(N, N, 3*N); 
Dyy = spalloc(N, N, 3*N); 

% For each interior point - compute its row in A
for k = K(:).'
    [i, j] = ind2sub(sz, k);
    if sz(1) > 1 % apply X derivative
        Dxx(k, sub2ind(sz, i+[1;0;-1], j+[0;0;0])) = ...
            ((C(i+1, j) + C(i, j)) * [1;-1; 0] / (X(i+1, j) - X(i, j)) - ...
             (C(i, j) + C(i-1, j)) * [0; 1;-1] / (X(i, j) - X(i-1, j))) / ...
             (X(i+1, j) - X(i-1, j));
    end

    if sz(2) > 1 % apply Y derivative
        Dyy(k, sub2ind(sz, i+[0;0;0], j+[1;0;-1])) = ...
            ((C(i, j+1) + C(i, j)) * [1;-1; 0] / (Y(i, j+1) - Y(i, j)) - ...
             (C(i, j) + C(i, j-1)) * [0; 1;-1] / (Y(i, j) - Y(i, j-1))) / ...
             (Y(i, j+1) - Y(i, j-1));
    end
end
A = Dxx + Dyy;
