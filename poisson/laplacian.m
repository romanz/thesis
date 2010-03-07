function [A, F, interior] = laplacian(sz, X, Y, C, F)
%% Laplacian discretization using sparse matrix
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
Ip = repmat(I, [1 3]);
Jp = repmat(J, [1 3]);
Kp = repmat(K, [1 3]);

% Build Laplacian sparse matrix:

if sz(1) > 1 % for X
    Kr = ind(I+1, J); % Left
    Kl = ind(I-1, J); % Right
    Dxx = ... % Laplacian stencil in X direction
     col(( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
     col(( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
    P = ind(Ip + Dp, Jp);
    Dxx = sparse(Kp, P, Dxx, N, N);
    % The denumerator is separated, for A to be symmetric
    P = ind(I, J);
    Mxx = sparse(P, P, (X(Kr) - X(Kl)), N, N);
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
    P = ind(Ip, Jp + Dp);
    Dyy = sparse(Kp, P, Dyy, N, N);
    % The denumerator is separated, for A to be symmetric
    P = ind(I, J);
    Myy = sparse(P, P, (Y(Ku) - Y(Kd)), N, N);
else
    Dyy = sparse(N, N);
    Myy = speye(N, N);
end
% Laplacian for Y direction is (pinv(Myy) * Dyy)

% pinv(Mxx) * Dxx + pinv(Myy) * Dyy = F
A = Myy * Dxx + Mxx * Dyy;
F = Mxx * Myy * F(:);
