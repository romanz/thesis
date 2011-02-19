function spheric
    U = @(X, Y) X.^2;
    Lu = @(X, Y) 6 .* X.^2 .* sin(Y);
    m = 4;
    x = linspace(log10(1), log10(2), 1+2^m);
    y = linspace(0, pi, 1+2^m);
    [sz, Xc, Yc, Xx, Yx, Xy, Yy, Xi, Yi] = create_grids(x, y);
    Ic = interior(sz + 2);    
    [Ju Iu Jd Id Jl Il Jr Ir] = boundary(Ic);
    K  = [Ju; Jd; Jl; Jr]; % (Xc, Yc)
    M = 0*K + 1;
    u = Xc(K) .* cos( Yc(K) );
    u = U(Xc(K), Yc(K));
    L = laplacian(Ic, Xc, Yc);    
    Cx = Xc.^2 .* sin(Yc);
    Cy = sin(Yc);
    L1 = L(Cx, Cy);
    f = shave(Lu(Xc, Yc), 1, 1);
    [A, f] = subst(L1, f, K, M, u);
    Phi = A \ f;
    Phi = reshape(Phi, sz);
    e = Phi - shave(U(Xc, Yc), 1, 1);
    norm(e(:), inf)
    figure(1); mesh(average(y, [1 1]/2), average(x, [1 1]/2), Phi)
    figure(2); mesh(average(y, [1 1]/2), average(x, [1 1]/2), e);
    save results
end

function A = shave(A, rows, cols)
    A = A(1+rows : end-rows, 1+cols : end-cols);
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
    function A = op(Cx, Cy)

        Dxx = ... % Laplacian stencil in X direction
         col(( Cx(Kr) + Cx(K) ) ./ (( X(Kr) - X(K) ))) * [1 -1 0] -  ...
         col(( Cx(K) + Cx(Kl) ) ./ (( X(K) - X(Kl) ))) * [0 1 -1];
        Dxx = sparse(Kp, Pi, Dxx, numel(K), N);
        
        Dyy = ... % Laplacian stencil in Y direction
         col(( Cy(Ku) + Cy(K) ) ./ (( Y(Ku) - Y(K) ))) * [1 -1 0] -  ...
         col(( Cy(K) + Cy(Kd) ) ./ (( Y(K) - Y(Kd) ))) * [0 1 -1];
        Dyy = sparse(Kp, Pj, Dyy, numel(K), N);
        
        L = Myy * Dxx + Mxx * Dyy;
        A = dinvM * L;
    end
    operator = @op;
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

function I = interior(sz)
    % Array of true values
    I = true(sz);
    % Exclude boundary points' indices (special handling of 1D)
    I([1 sz(1)], :) = false;
    I(:, [1 sz(2)]) = false; 
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
