function [A] = advection(interior, X, Y, Vx, Vy)
%% Advection gradient discretization using sparse matrix.
sz = size(interior);
N = numel(interior);

ind = @(I, J) sub2ind(sz, I, J);
K = find(interior); % Fill only interior points
K = K(:);
[I, J] = ind2sub(sz, K);
Dp = repmat([1 0 -1], [numel(K) 1]); 
% NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
Ip = repmat(I, [1 3]);
Jp = repmat(J, [1 3]);
Kp = repmat((1:numel(K))', [1 3]); % interior variables' indices, for 1D stencil

% for X
Vx = Vx(:, 2:end-1); % remove velocity at ghost-points
Kr = ind(I+1, J); % Left from K
Kl = ind(I-1, J); % Right from K
P = ind(Ip + Dp, Jp); % column indices
Vx = average(Vx, [1; 1]/2);
Dx = spdiag( Vx(:) ./ (X(Kr) - X(Kl)) ) * ...
     sparse(Kp, P, ones(size(K)) * [1 0 -1], numel(K), N);

% for Y
Vy = Vy(2:end-1, :); % remove velocity at ghost-points
Ku = ind(I, J+1); % Up from K
Kd = ind(I, J-1); % Down from K
P = ind(Ip, Jp + Dp); % column indices
Vy = average(Vy, [1, 1]/2);
Dy = spdiag( Vy(:) ./ ((Y(Ku) - Y(Kd)) .* X(K)) ) * ...
     sparse(Kp, P, ones(size(K)) * [1 0 -1], numel(K), N);        

A = Dx + Dy;
