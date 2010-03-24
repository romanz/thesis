function [VG, interior] = gradient(interior, X, Y, Vx, Vy, method)
%% Gradient discretization using sparse matrix
sz = size(interior);

ind = @(I, J) sub2ind(sz, I, J);
K = find(interior); % Fill only interior points
K = K(:);
[I, J] = ind2sub(sz, K);
Dp = repmat([1 0 -1], [numel(K) 1]); 
% NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
Ip = repmat(I, [1 3]);
Jp = repmat(J, [1 3]);
Kp = repmat(K, [1 3]); % interior variables' indices, for 1D stencil

switch method
    case 'upstream'
        Dx = sparse(N, N);
        if sz(1) > 1 % for X
            Kr = ind(I+1, J); % Left
            Kl = ind(I-1, J); % Right
            P = ind(Ip + Dp, Jp); % column indices
            Dx_r = sparse(K, K, 1./(X(Kr) - X(K)), N, N) * ...
                   sparse(Kp, P, ones(size(K)) * [1 -1 0], N, N);
            Dx_l = sparse(K, K, 1./(X(K) - X(Kl)), N, N) * ...
                   sparse(Kp, P, ones(size(K)) * [0 1 -1], N, N);
            % Find upstream direction
            Sx = sign(convn(Vx, [1; 1]/2, 'valid'));
            Dx = sparse(K, K, (Sx > 0) .* Vx(1:end-1, :), N, N) * Dx_l + ...
                 sparse(K, K, (Sx < 0) .* Vx(2:end,   :), N, N) * Dx_r;
        end
        Dy = sparse(N, N);
        if sz(2) > 1 % for Y
            Ku = ind(I, J+1); % Up
            Kd = ind(I, J-1); % Down
            P = ind(Ip, Jp + Dp); % column indices
            Dy_u = sparse(K, K, 1./(Y(Ku) - Y(K)), N, N) * ...
                 sparse(Kp, P, ones(size(K)) * [1 -1 0], N, N);
            Dy_d = sparse(K, K, 1./(Y(K) - Y(Kd)), N, N) * ...
                 sparse(Kp, P, ones(size(K)) * [0 1 -1], N, N);
            % Find upstream direction
            Sy = sign(convn(Vy, [1, 1]/2, 'valid'));
            Dy = sparse(K, K, (Sy > 0) .* Vy(:, 1:end-1), N, N) * Dy_d + ...
                 sparse(K, K, (Sy < 0) .* Vy(:, 2:end  ), N, N) * Dy_u;
        end
        VG = Dx + Dy;
    case 'centered'
        Dx = sparse(N, N);
        if sz(1) > 1 % for X
            Kr = ind(I+1, J); % Left
            Kl = ind(I-1, J); % Right
            P = ind(Ip + Dp, Jp); % column indices
            Vx = convn(Vx, [1; 1]/2, 'valid');
            Dx = sparse(K, K, Vx(:)./(X(Kr) - X(Kl)), N, N) * ...
                 sparse(Kp, P, ones(size(K)) * [1 0 -1], N, N);
        end
        Dy = sparse(N, N);
        if sz(2) > 1 % for Y
            Ku = ind(I, J+1); % Up
            Kd = ind(I, J-1); % Down
            P = ind(Ip, Jp + Dp); % column indices
            Vy = convn(Vy, [1, 1]/2, 'valid');
            Dy = sparse(K, K, Vy(:)./(Y(Ku) - Y(Kd)), N, N) * ...
                 sparse(Kp, P, ones(size(K)) * [1 0 -1], N, N);
        end
        VG = Dx + Dy;
        
    otherwise
        error('Unsupported %s!', method)
end
