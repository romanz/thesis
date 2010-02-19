function [A, I] = laplacian(X, Y, sz)
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
A = spalloc(N, N, 5*N); % Preallocate sparse matrix.
P = [0 -1 1; 1 0 -1; -1 1 0]; % Used for interpolation
% For each interior point - compute its row in A
for k = K(:).'
    if sz(1) > 1
        [i, j] = ind2sub(sz, k);
        i = i + [-1; 0; 1];
        j = j + [ 0; 0; 0];    
        m = sub2ind(sz, i, j);
        x = X(m);
        x = P * x; % [x3-x2; x1-x3; x2-x1]
        x = x / prod(x);
        A(k, m) = A(k, m) + x.';
    end

    if sz(2) > 1
        [i, j] = ind2sub(sz, k);
        i = i + [ 0; 0; 0];
        j = j + [-1; 0; 1];
        m = sub2ind(sz, i, j);
        y = Y(m);
        y = P * y; % [y3-y2; y1-y3; y2-y1]
        y = y / prod(y);
        A(k, m) = A(k, m) + y.';
    end
end
A = -2*A; % Fix missing factor and sign-reverse.
