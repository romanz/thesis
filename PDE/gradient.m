function [A, interior] = gradient(sz, X, Y)
%% Gradient discretization using sparse matrix
N = prod(sz);
interior = true(sz);

% Compute interior points' indices (special handling of 1D)
if sz(1) > 1,   interior([1 sz(1)], :) = false; end
if sz(2) > 1,   interior(:, [1 sz(2)]) = false; end

ind = @(I, J) sub2ind(sz, I, J);
K = find(interior); % Fill only interior points
K = K(:);
[I, J] = ind2sub(sz, K);
Dp = repmat([1 0 -1], [numel(K) 1]); 
% NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
Ip = repmat(I, [1 3]);
Jp = repmat(J, [1 3]);
Kp = repmat(K, [1 3]); % interior variables' indices, for 1D stencil

% Build 2D sparse matrix A:

if sz(1) > 1 % for X
    Kr = ind(I+1, J); % Left
    Kl = ind(I-1, J); % Right
    Gx = ones(size(K)) * [1 0 -1]; 
    P = ind(Ip + Dp, Jp); % column indices
    Gx = sparse(Kp, P, Gx, N, N);
    Mx = sparse(K, K, (X(Kr) - X(Kl)), N, N);
else
    Gx = sparse(N, N);
    Mx = speye(N, N);
end

if sz(2) > 1 % for Y
    Ku = ind(I, J+1); % Up
    Kd = ind(I, J-1); % Down
    Gy = ones(size(K)) * [1 0 -1];
    P = ind(Ip, Jp + Dp); % column indices
    Gy = sparse(Kp, P, Gy, N, N);
    My = sparse(K, K, (Y(Ku) - Y(Kd)), N, N);
else
    Gy = sparse(N, N);
    My = speye(N, N);
end

A = dinv(Mx) * Gx + dinv(My) * Gy;
