function A = laplacian(interior, X, Y, C)
%% Laplacian discretization using sparse matrix
sz = size(interior);
N = prod(sz);
if nargin < 4
    C = ones(sz); % Assume uniform C
end

% Fill only interior points
M = nnz(interior);
[I, J] = ind2sub(sz, find(interior));
Dp = repmat([1 0 -1], [M 1]); 
% NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
Ip = repmat(I, [1 3]);
Jp = repmat(J, [1 3]);
Kp = repmat((1:M)', [1 3]); % interior variables' indices, for 1D Laplacian stencil

% Build 2D sparse matrix L = "Dxx/Mxx" + "Dyy/Myy":
K = interior;
Kr = shift(K, [+1 0]); % Left
Kl = shift(K, [-1 0]); % Right
Ku = shift(K, [0 +1]); % Up
Kd = shift(K, [0 -1]); % Down
ind = @(I, J) sub2ind(sz, I, J);

% for X
    Dxx = ... % Laplacian stencil in X direction
     (( (X(Kr) + X(K))/2 ).^2 .* ( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
     (( (X(K) + X(Kl))/2 ).^2 .* ( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
    P = ind(Ip + Dp, Jp); % column indices
    Dxx = sparse(Kp, P, Dxx, M, N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Mxx = sparse(1:M, 1:M, X(K).^2 .* (X(Kr) - X(Kl)), M, M);
% Laplacian for X direction is (pinv(Mxx) * Dxx)

% for Y
    Dyy = ... % Laplacian stencil in Y direction
     (sin( (Y(Ku) + Y(K))/2 ) .* ( C(Ku) + C(K) ) ./ (( Y(Ku) - Y(K) ))) * [1 -1 0] -  ...
     (sin( (Y(K) + Y(Kd))/2 ) .* ( C(K) + C(Kd) ) ./ (( Y(K) - Y(Kd) ))) * [0 1 -1];
    P = ind(Ip, Jp + Dp); % column indices
    Dyy = sparse(Kp, P, Dyy, M, N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Myy = sparse(1:M, 1:M, X(K).^2 .* sin(Y(K)) .* (Y(Ku) - Y(Kd)), M, M);
% Laplacian for Y direction is (pinv(Myy) * Dyy)

% % The actual Laplacian matrix is:
% A = dinv(Mxx) * Dxx + dinv(Myy) * Dyy;

% Pre-multiply it by (Mxx * Myy) for symmetry of A.
M = Mxx * Myy;
L = Myy * Dxx + Mxx * Dyy;
A = dinv(M) * L;
return

dinv(Mxx) * Dxx
dinv(Myy) * Dyy
