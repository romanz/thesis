function [A] = advection_upwind(interior, X, Y, Vx, Vy, gridVx, gridVy)
%% Advection gradient upwind discretization using sparse matrix.
sz = size(interior);
N = numel(interior);

ind = @(I, J) sub2ind(sz, I, J);
K = find(interior); % Fill only interior points
K = K(:);
[I, J] = ind2sub(sz, K); % Interior subscript values

% for X
Sx = Vx .* (gridVx.X .^ 2); % Compute flux
Sx = Sx(:, 2:end-1); % remove velocity at ghost-points
Sx = average(Sx, [1; 1]) > 0; % 1 = positive flux, 0 = negative flux
Vx = Vx(1:end-1, 2:end-1) .* Sx + Vx(2:end, 2:end-1) .* ~Sx;
Kl = ind(I-Sx(:)  , J); % Left point index
Kr = ind(I-Sx(:)+1, J); % Right point index
vDx = spdiag( Vx(:) ./ (X(Kr) - X(Kl)) ) * ...
     (sparse(1:numel(K), Kr, 1, numel(K), N) - sparse(1:numel(K), Kl, 1, numel(K), N));

% for Y
Sy = Vy(2:end-1, :); % remove velocity at ghost-points
Sy = average(Sy, [1, 1]) > 0; % 1 = positive flux, 0 = negative flux
Vy = Vy(2:end-1, 1:end-1) .* Sy + Vy(2:end-1, 2:end) .* ~Sy;
Kd = ind(I, J-Sy(:));
Ku = ind(I, J-Sy(:)+1);
vDy = spdiag( Vy(:) ./ ((Y(Ku) - Y(Kd)) .* X(K)) ) * ...
    (sparse(1:numel(K), Ku, 1, numel(K), N) - sparse(1:numel(K), Kd, 1, numel(K), N));

A = vDx + vDy;
