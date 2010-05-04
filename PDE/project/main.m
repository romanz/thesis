% The main solver code.
function main
    clc;
    % Physical quantities
    alpha = 0; % Peclet number
    beta  = 0; % Applied electric field
    vel = 0; % Particle velocity 
    
    % Grid creation
    x = linspace(-1, 1, 5);
    y = linspace(-1, 1, 5);
    [sz, Xc, Yc, Xx, Yx, Xy, Yy, Xi, Yi] = create_grids(x, y, 1);
    Ic = interior(sz + 2);
    [Ju Iu Jd Id Jl Il Jr Ir] = boundary(Ic);
    
    % Guess initial values for all variables
    Phi = zeros(sz) + randn(sz)/100;
    C = ones(sz + 2);
    C(Ic) = 1 + randn(sz)/100;
    Vx = zeros(sz - [1 0]);
    Vy = zeros(sz - [0 1]);
    P = zeros(sz);
    
    % Iterations
    iters = [4 2];
    L = laplacian(Ic, Xc, Yc);
    for iter = 1:1000
        %%% Laplace equation (for Phi)
        A = laplacian(Ic, Xc, Yc, C);
        K = [Ju Iu; Jd Id; Jl Il; Jr Ir];
        M = [0*Ju+1 0*Iu-1; 0*Jd+1 0*Id-1; 0*Jl+1 0*Il; 0*Jr+1 0*Ir-1];
        u = [0*[Ju; Jd]; -log(col(C(Il))); 0*Jr];
        [A, f] = subst(A, zeros(sz), K, M, u);
        Phi = iterate(Phi, A, f, iters(1));
        %%% Diffusion-advection (for C)
        M = [0*Ju+1 0*Iu-1; 0*Jd+1 0*Id-1; 0*Jl+1 0*Il; 0*Jr+1 0*Ir];
        u = [0*[Ju; Jd]; exp(-col(Phi(1, 1:end))); 0*Jr+1];
        A = L;
        [A, f] = subst(A, zeros(sz), K, M, u);
        C(Ic) = iterate(C(Ic), A, f, iters(2));
        C(Jl) = exp(-col(Phi(1, 1:end)));        
        C(Jr) = 1;
        C(Ju) = C(Iu);
        C(Jd) = C(Id);
        
        %%% Stokes (for Vx, Vy and P)        
        % P = P - mean(P(:)); % re-normalize P
    end
    % Results
    Phi, C-1, Vx, Vy, P

    save results
end

function x = iterate(x, A, f, N)
    sz = size(x);
    x = x(:);
    for i = 1:N
        r = f - A * x;
        x = x + dinv(A) * r;
    end
    x = reshape(x, sz);
end

function [Ju Iu Jd Id Jl Il Jr Ir] = boundary(I)
    Jl = find(  shift(I, [-1, 0]) & ~I ); 
    Il = find( ~shift(I, [+1, 0]) &  I ); 
    Jr = find(  shift(I, [+1, 0]) & ~I ); 
    Ir = find( ~shift(I, [-1, 0]) &  I ); 
    Jd = find(  shift(I, [0, -1]) & ~I ); 
    Id = find( ~shift(I, [0, +1]) &  I ); 
    Ju = find(  shift(I, [0, +1]) & ~I ); 
    Iu = find( ~shift(I, [0, -1]) &  I ); 
end

function I = shift(I, Q)

    sz = size(I);
    Qp = max(Q, 0); % positive shifts
    Qn = min(Q, 0); % negative shifts

    I = [zeroes(Q(1), size(I, 2)); ...
         I(1-Qn(1):end-Qp(1), :); ...
         zeroes(-Q(1), size(I, 2))];

    I = [zeroes(size(I, 1), Q(2)), ...
         I(:, 1-Qn(2):end-Qp(2)), ...
         zeroes(size(I, 1), -Q(2))];

    function Z = zeroes(n, m)
        Z = false(min(max([n m], 0), sz));
    end

end

% Substitude Mx(K)=u into Ax=f, removing eliminated variables.
% The first coordinate for each M's row is eliminated.
% In addition, all zero columns are eliminated from A.
function [A, f] = subst(A, f, K, M, u)

    P = K(:, 1); % P is to-be-eliminated variable vector
    K(:, 1) = []; % K now contains dependent variables
    S = diag(1 ./ M(:, 1)); % Rescaling according to first coordinate

    u = S * u(:);
    f = f(:) - A(:, P) * u;

    n = size(A, 2); % n = # of variables.
    
    M = S * M(:, 2:end);
    T = sparse(repmat(P, 1, size(K, 2)), K, M, n, n);
    A = A - A * T;
    A(:, P) = [];

    A(:, ~any(A, 1)) = [];
end

function I = interior(sz)
    % Array of true values
    I = true(sz);
    % Exclude boundary points' indices (special handling of 1D)
    I([1 sz(1)], :) = false;
    I(:, [1 sz(2)]) = false; 
end

% Laplacian discretization using sparse matrix
function A = laplacian(interior, X, Y, C)    
    sz = size(interior);
    N = prod(sz);
    if nargin < 4
        C = ones(sz); % Assume uniform C
    end

    ind = @(I, J) sub2ind(sz, I, J);
    K = find(interior); % Fill only interior points
    K = K(:);
    [I, J] = ind2sub(sz, K);
    Dp = repmat([1 0 -1], [numel(K) 1]); 
    % NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
    Ip = repmat(I, [1 3]);
    Jp = repmat(J, [1 3]);

    % Build 2D sparse matrix L = "Dxx/Mxx" + "Dyy/Myy":

    % for X
    Kr = ind(I+1, J); % Left
    Kl = ind(I-1, J); % Right
    Dxx = ... % Laplacian stencil in X direction
     col(( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
     col(( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
    P = ind(Ip + Dp, Jp); % column indices
    Dxx = sparse(repmat((1:numel(K))', [1 3]), P, Dxx, numel(K), N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Mxx = sparse(1:numel(K), 1:numel(K), (X(Kr) - X(Kl)), numel(K), numel(K));
    % Laplacian for X direction is (pinv(Mxx) * Dxx)

    % for Y
    Ku = ind(I, J+1); % Up
    Kd = ind(I, J-1); % Down
    Dyy = ... % Laplacian stencil in Y direction
     col(( C(Ku) + C(K) ) ./ (( Y(Ku) - Y(K) ))) * [1 -1 0] -  ...
     col(( C(K) + C(Kd) ) ./ (( Y(K) - Y(Kd) ))) * [0 1 -1];
    P = ind(Ip, Jp + Dp); % column indices
    Dyy = sparse(repmat((1:numel(K))', [1 3]), P, Dyy, numel(K), N);
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Myy = sparse(1:numel(K), 1:numel(K), (Y(Ku) - Y(Kd)), numel(K), numel(K));
    % Laplacian for Y direction is (pinv(Myy) * Dyy)

    % % The actual Laplacian matrix is:
    % A = dinv(Mxx) * Dxx + dinv(Myy) * Dyy;

    % Pre-multiply it by (Mxx * Myy) for symmetry of A.
    M = Mxx * Myy;
    L = Myy * Dxx + Mxx * Dyy;
    A = dinv(M) * L;
end

% Diffusion advection discretization using sparse matrix
function [A, interior] = advection(interior, X, Y, Vx, Vy, method)    
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
    Kp = repmat(K, [1 3]); % interior variables' indices, for 1D stencil

    switch method
        case 'upwind'
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

        case 'central'
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

        otherwise
            error('Unsupported %s!', method)
    end

    A = Dx + Dy;
end

function [sz, Xc, Yc, Xx, Yx, Xy, Yy, Xi, Yi] = create_grids(x, y, verbose)
    x = x(:); y = y(:); sz = [numel(x) numel(y)] - 1;
    % extended grid (for ghost points)
    xg = [2*x(1) - x(2); x; 2*x(end) - x(end-1)];
    yg = [2*y(1) - y(2); y; 2*y(end) - y(end-1)];
    % cell-centered coordinates
    xc = average(xg, [1; 1]/2);
    yc = average(yg, [1; 1]/2);
    % NDGRID convention
    [Xc, Yc] = ndgrid(xc, yc);
    [Xx, Yx] = ndgrid(xg(2:end-1), yc);
    [Xy, Yy] = ndgrid(xc, yg(2:end-1));
    [Xi, Yi] = ndgrid(xc(2:end-1), yc(2:end-1));

    if nargin > 2 && verbose        
        plot(Xc(:), Yc(:), '.', Xi(:), Yi(:), 'o', ...
             Xx(:), Yx(:), '>', Xy(:), Yy(:), '^');
        set(gca, 'XTick', x, 'YTick', y, 'XGrid', 'on', 'YGrid', 'on');
        axis square; 
    end
end

function A = average(A, h)
    A = convn(A, h, 'valid');
end

% R = pseudo-inverse of diag(A).
function R = dinv(A)
    d = full(diag(A));
    I = find(d);
    d(I) = 1 ./ d(I);
    N = numel(d);
    R = sparse(1:N, 1:N, d);
end

% Converts x to a column vector = x(:)
function x = col(x)
    x = x(:);
end
