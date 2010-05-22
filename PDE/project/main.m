% The main solver code.
function main
    % clc;
    % Physical quantities
    alpha = 1; % Peclet number
    beta = 0.1; % Applied electric field
    vel = 0.1; % Particle velocity 
    
    % Grid creation
    x = logspace(log10(1), log10(2), 31);
    y = linspace(0, 1, 21);
    [sz, Xc, Yc, Xx, Yx, Xy, Yy, Xi, Yi] = create_grids(x, y);
    Ic = interior(sz + 2);
    [Ju Iu Jd Id Jl Il Jr Ir] = boundary(Ic);
    K  = [Ju Iu; Jd Id; Jl Il; Jr Ir]; % (Xc, Yc)
    
    [Jux Iux Jdx Idx Jlx Ilx Jrx Irx] = boundary(interior(sz + [1 2]));
    Kx = [Jux Iux; Jdx Idx; Jlx Ilx; Jrx Irx]; % (Xx, Yx)
    Mx = [0*Jux+1, 0*Iux-1; 0*Jdx+1, 0*Idx-1; ...
          0*Jlx+1, 0*Ilx; 0*Jrx+1 0*Irx];
    % Up/Down: Neumann = 0/0.
    % Left/Right: Dirichlet = 0/cos(theta)
    
    [Juy Iuy Jdy Idy Jly Ily Jry Iry] = boundary(interior(sz + [2 1]));
    Ky = [Juy Iuy; Jdy Idy; Jly Ily; Jry Iry]; % (Xy, Yy)
    My = [0*Juy+1, 0*Iuy; 0*Jdy+1, 0*Idy; ...
          0*[Jly, Ily; Jry, Iry]+0.5];
    % Up/Down: Dirichlet = 0/0
    % Left/Right: Average Dirichlet = Slip/sin(theta)
    
    % Guess initial values for all variables
    Phi = randn(sz);
    C = rand(sz);
    Vx = randn(sz - [1 0]);
    Vy = randn(sz - [0 1]);
    P = randn(sz);
    
    % Iterations
    iters = [1 2 10];
    L = laplacian(Ic, Xc, Yc);
    L1 = L();
    [Gx_p, L_vx0, Gx_vx0] = stokes1(1, Xx, Yx, Xi, Yi);
    [Gy_p, L_vy0, Gy_vy0] = stokes1(2, Xy, Yy, Xi, Yi);
    
    [L_vx] = subst_lhs(L_vx0, Kx, Mx);
    [L_vy] = subst_lhs(L_vy0, Ky, My);
    [Gx_vx] = subst_lhs(Gx_vx0, Kx, Mx);
    [Gy_vy] = subst_lhs(Gy_vy0, Ky, My);

    A1 = [L_vx, sparse(size(L_vx, 1), size(L_vy, 2)), -Gx_p; ...
         sparse(size(L_vy, 1), size(L_vx, 2)), L_vy, -Gy_p; ...
         Gx_vx, Gy_vy, sparse(prod(sz), prod(sz))];
    [Mr, Mb] = stokes_vanka_redblack(sz, A1);

    VxL = zeros(1, size(Vx, 2));
    VxR = zeros(1, size(Vx, 2)) - vel*cos(Yx(end, 2:end-1)*pi);
    VyU = zeros(size(Vy, 1), 1);
    VyD = zeros(size(Vy, 1), 1);

    Cl = exp(-Phi(1, :));
    Cr = 0*Cl + exp(-0);
    
    for iter = 1:10000
        %%% Laplace equation (for Phi)
        if iters(3)
            C0 = [Cl; C; Cr]; % Expand C with ghost points
            A = L([C0(:, 1), C0, C0(:, end)]);
            % Up/Down: symmetry - Neumann. Left: Dirichlet. Right: Neumann (field).
            M = [0*Ju+1 0*Iu-1; 0*Jd+1 0*Id-1; 0*Jl+1 0*Il; 0*Jr+1 0*Ir-1];
            u = [0*[Ju; Jd]; -log(col(C(1, :))); ...
                beta * (Xc(Jr) - Xc(Ir)) .* cos(Yc(Jr) * pi)];
            [A, f] = subst(A, zeros(sz), K, M, u);
            [Phi, res3] = iterate(Phi, A, f, iters(3));
            Cl = exp(-Phi(1, :));
        end
        %%% Diffusion-advection (for C)        
        if iters(2)
            % Up/Down: symmetry - Neumann. Left: Dirichlet. Right: Dirichlet.
            M = [0*Ju+1 0*Iu-1; 0*Jd+1 0*Id-1; 0*Jl+1 0*Il; 0*Jr+1 0*Ir];
            u = [0*[Ju; Jd]; Cl(:); Cr(:)];

            A = advection(Ic, Xc, Yc, [VxL; Vx; VxR], [VyD, Vy, VyU], 'central');        
            [A, f] = subst(L1 - alpha*A, zeros(sz), K, M, u);
            [C, res2] = iterate(C, A, f, iters(2));
        end        
        %%% Stokes equation (for Vx, Vy, P)
        if iters(1)
            Fx = shave(Xx, 1, 1) * 0;
            Fy = shave(Yy, 1, 1) * 0;
            div = zeros(sz);

            %    Neumann (U/D)     Dirichlet (L/R)
            ux = [[0*Jux; 0*Jdx]; [VxL(:); VxR(:)]];
            [Fx] = subst_rhs(L_vx0, Fx, Kx, Mx, ux);
            [div] = subst_rhs(Gx_vx0, div, Kx, Mx, ux);
            
            %    Dirichlet (U/D)   Av. Dirichlet (L/R)  [TODO]
            uy = [[VyU(:); VyD(:)]; [0*Jly; 0*Jry]];
            [Fy] = subst_rhs(L_vy0, Fy, Ky, My, uy);
            [div] = subst_rhs(Gy_vy0, div, Ky, My, uy);

            f1 = [Fx; Fy; div];

            z = [Vx(:); Vy(:); P(:)];
            for k = 1:iters(1)
                r = f1 - A1*z;    z = z + Mr*r; 
                r = f1 - A1*z;    z = z + Mb*r;
            end
            res1 = r;
            Vx = reshape(z(1:numel(Fx)), sz - [1 0]);
            Vy = reshape(z((1:numel(Fy)) + numel(Fx)), sz - [0 1]);
            P = reshape(z((1 + numel(Fx) + numel(Fy)):end), sz);
            P = P - mean(P(:)); % re-normalize P
        end
    end
    
    % Results
%     Phi, C, Vx, Vy, P
    figure(1); 
    mesh(Phi); colorbar; title('\Phi')
    
    figure(2); 
    mesh(C); colorbar; title('C')

    figure(3); 
    quiver(Xi, Yi, ...
        average([VxL; Vx; VxR], [1;1]/2), ...
        average([VyD, Vy, VyU], [1 1]/2), 0); title('V')
    axis([1 2 0 1])

    figure(4); 
    imagesc(average(x, [1 1]/2), average(y, [1 1]/2), P.'); 
    colorbar; title('P')
    set(gca, 'YDir', 'normal');
    save results
    norm(res1, inf)
    norm(res2, inf)
    norm(res3, inf)
end

function A = shave(A, rows, cols)
    A = A(1+rows : end-rows, 1+cols : end-cols);
end

function [x, r] = iterate(x, A, f, N)
    sz = size(x);
    x = x(:);
    J = dinv(A);
    for i = 1:N
        r = f - A * x;
        x = x + J * r;
    end
    x = reshape(x, sz);
end

function [M1, M2] = stokes_vanka_redblack(sz, A)
    K = 1:prod(sz);
    K = K(:);
    [I, J] = ind2sub(sz, K);
    
    % Same indices for variables and equations.
    % CR: make this offset hack better.
    V = [[index(sz - [1 0], I-1, J), index(sz - [1 0], I, J)], ...
         [index(sz - [0 1], I, J-1), index(sz - [0 1], I, J)] + prod(sz - [1 0]), ...
         K + prod(sz - [1 0]) + prod(sz - [0 1])];
    
    K = logical(mod(I - J, 2)); 
    V1 = arr2cell(V( K, :)); % odd
    V2 = arr2cell(V(~K, :)); % even
    
    M1 = vanka(A, V1, V1);
    M2 = vanka(A, V2, V2);

    function C = arr2cell(A)
        C = cell(size(A, 1), 1);
        for k = 1:numel(C)
            a = A(k, :);
            a = a(~isnan(a));
            C{k} = a(:);
        end
    end

    function K = index(sz, I, J)
        Q = (1 <= I) & (I <= sz(1)) & (1 <= J) & (J <= sz(2));
        K = nan(size(Q));
        K(Q) = sub2ind(sz, I(Q), J(Q));
    end

end

% Vanka-type smoother construction.
function [M] = vanka(A, V, E)

    % V is a cell array, whose k-th entry contains k-th subdomain indices
    nonzeros = 0;
    for k = 1:numel(V)
        nonzeros = nonzeros + numel(V{k}).^2;
    end
    % Preallocate memory for non-zeroes
    inverses = zeros(nonzeros, 1);
    indices = zeros(size(inverses));

    sz = size(A); % Pre-compute A size
    offset = 0;
    for k = 1:numel(V)
         % Process current index set:
        J = col(V{k});
        I = col(E{k});
        N = numel(J); 
        % Invert current submatrix
        inverses( offset + (1:N^2) ) = inv( A(I, J) );
        % Save its indices at A
        indices( offset + (1:N^2) ) = index_ndgrid(sz, J, I);
        offset = offset + N^2;
    end
    % Construct sparse matrix M for pre-conditioning
    M = mksparse(sz, indices, inverses);

    % r = f - Ax
    % x' = x + M r 
    % x' = (I - MA)x + Mf
    % x' = Tx + d
end

function K = index_ndgrid(sz, I, J)
    [I, J] = ndgrid(I, J);
    K = sub2ind(sz, I, J);
end
function S = mksparse(sz, indices, values)
    [I, J] = ind2sub(sz, indices);
    S = sparse(I, J, values, sz(1), sz(2));
end

% Construct Stokes equation for specified dimension.
function [Gp, Lv, Gv] = stokes1(dim, Xs, Ys, X, Y)
    sz = size(0*Xs + 0*Ys);
    I = interior(sz);
    Lv = laplacian(I, Xs, Ys);
    Lv = Lv();
    % full(Lv)
    dir = (1:numel(sz)) == dim;
    K = I | shift(I, -dir); % add first row/column
    Gv = gradient(K, Xs, Ys, dir);
    % full(Gv)
    I = true(size(0*X+0*Y));
    K = I & shift(I, -dir); % remove last row/column
    Gp = gradient(K, X, Y, dir);
    % full(Gp)
end

% Sparse gradient operator.
function G = gradient(K, X, Y, dir)
    sz = size(K);
    [I, J] = ind2sub(sz, find(K));
    K1 = sub2ind(sz, I, J);
    K2 = sub2ind(sz, I + dir(1), J + dir(2));
    K0 = sub2ind(sz - dir, I, J);
    switch find(dir)
        case 1, D = X(K2) - X(K1);
        case 2, D = Y(K2) - Y(K1);
    end
    G = sparse([1:numel(K0), 1:numel(K0)], [K1 K2], 1./[-D D], ...
        numel(K0), prod(sz));
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
    f = subst_rhs(A, f, K, M, u);
    A = subst_lhs(A, K, M);
end

function [A] = subst_lhs(A, K, M)

    P = K(:, 1); % P is to-be-eliminated variable vector
    K(:, 1) = []; % K now contains dependent variables
    S = diag(1 ./ M(:, 1)); % Rescaling according to first coordinate

    n = size(A, 2); % n = # of variables.
    
    M = S * M(:, 2:end);
    T = sparse(repmat(P, 1, size(K, 2)), K, M, n, n);
    A = A - A * T;
    A(:, P) = [];

    A(:, ~any(A, 1)) = [];
end

function [f] = subst_rhs(A, f, K, M, u)

    P = K(:, 1); % P is to-be-eliminated variable vector
    S = diag(1 ./ M(:, 1)); % Rescaling according to first coordinate

    u = S * u(:);
    f = f(:) - A(:, P) * u;
end

function I = interior(sz)
    % Array of true values
    I = true(sz);
    % Exclude boundary points' indices (special handling of 1D)
    I([1 sz(1)], :) = false;
    I(:, [1 sz(2)]) = false; 
end

% Laplacian discretization using sparse matrix
function operator = laplacian(interior, X, Y)    
    sz = size(interior);
    N = prod(sz);

    ind = @(I, J) sub2ind(sz, I, J);
    K = find(interior); % Fill only interior points
    K = K(:);
    [I, J] = ind2sub(sz, K);
    Dp = repmat([1 0 -1], [numel(K) 1]); 
    % NOTE: The stencil is constructed as: [forward, middle, backward] coefficients.
    Ip = repmat(I, [1 3]);
    Jp = repmat(J, [1 3]);
    Kp = repmat((1:numel(K))', [1 3]);
    % Build 2D sparse matrix L = "Dxx/Mxx" + "Dyy/Myy":

    % for X
    Kr = ind(I+1, J); % Left
    Kl = ind(I-1, J); % Right
    Pi = ind(Ip + Dp, Jp); % column indices
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Mxx = sparse(1:numel(K), 1:numel(K), (X(Kr) - X(Kl)), numel(K), numel(K));
    % Laplacian for X direction is (pinv(Mxx) * Dxx)

    % for Y
    Ku = ind(I, J+1); % Up
    Kd = ind(I, J-1); % Down
    Pj = ind(Ip, Jp + Dp); % column indices
    % The denumerator is separated, for A to be symmetric
    % (since the original operator is self-adjoint).
    Myy = sparse(1:numel(K), 1:numel(K), (Y(Ku) - Y(Kd)), numel(K), numel(K));
    % Laplacian for Y direction is (pinv(Myy) * Dyy)

    % % The actual Laplacian matrix is:
    % A = dinv(Mxx) * Dxx + dinv(Myy) * Dyy;

    % Pre-multiply it by (Mxx * Myy) for symmetry of A.
    M = Mxx * Myy;
    dinvM = dinv(M);
    function A = op(C)
        if nargin < 1
            C = ones(sz); % Assume uniform C
        end

        Dxx = ... % Laplacian stencil in X direction
         col(( C(Kr) + C(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
         col(( C(K) + C(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
        Dxx = sparse(Kp, Pi, Dxx, numel(K), N);
        
        Dyy = ... % Laplacian stencil in Y direction
         col(( C(Ku) + C(K) ) ./ (( Y(Ku) - Y(K) ))) * [1 -1 0] -  ...
         col(( C(K) + C(Kd) ) ./ (( Y(K) - Y(Kd) ))) * [0 1 -1];
        Dyy = sparse(Kp, Pj, Dyy, numel(K), N);
        
        L = Myy * Dxx + Mxx * Dyy;
        A = dinvM * L;
    end
    operator = @op;
end

% Diffusion advection discretization using sparse matrix
function A = advection(interior, X, Y, Vx, Vy, method)    
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
    A = A(interior, :);
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
