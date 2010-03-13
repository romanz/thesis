function [L, M, interior] = laplacian(sz, X, Y, C)
%% Laplacian discretization using sparse matrix
N = prod(sz);
interior = true(sz);

if nargin < 4
    C = ones(sz); % Assume uniform C
end

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
Kp = repmat(K, [1 3]); % interior variables' indices, for 1D Laplacian stencil

% Build 2D sparse matrix L = "Dxx/Mxx" + "Dyy/Myy":

if sz(1) > 1 % for X
    Kr = ind(I+1, J); % Left
    Kl = ind(I-1, J); % Right
    Dxx = ... % Laplacian stencil in X direction
     col(( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
     col(( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
    P = ind(Ip + Dp, Jp); % column indices
    Dxx = sparse(Kp, P, Dxx, N, N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Mxx = sparse(K, K, (X(Kr) - X(Kl)), N, N);
else
    Dxx = sparse(N, N);
    Mxx = speye(N, N);
end
% Laplacian for X direction is (pinv(Mxx) * Dxx)

if sz(2) > 1 % for Y
    Ku = ind(I, J+1); % Up
    Kd = ind(I, J-1); % Down
    Dyy = ... % Laplacian stencil in Y direction
     col(( C(Ku) + C(K) ) ./ (( Y(Ku) - Y(K) ))) * [1 -1 0] -  ...
     col(( C(K) + C(Kd) ) ./ (( Y(K) - Y(Kd) ))) * [0 1 -1];
    P = ind(Ip, Jp + Dp); % column indices
    Dyy = sparse(Kp, P, Dyy, N, N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Myy = sparse(K, K, (Y(Ku) - Y(Kd)), N, N);
else
    Dyy = sparse(N, N);
    Myy = speye(N, N);
end
% Laplacian for Y direction is (pinv(Myy) * Dyy)

% % The actual Laplacian matrix is:
% A = dinv(Mxx) * Dxx + dinv(Myy) * Dyy;

% Pre-multiply it by (Mxx * Myy) for symmetry of A.
M = Mxx * Myy;
L = Myy * Dxx + Mxx * Dyy;
