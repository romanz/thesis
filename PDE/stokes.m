% 2D Stokes equation discretization for specific dimension
% {X,Y} are cell centers coordinates (incl. the boundary cells)
function [A, f] = stokes(sz, X, Y, Fx, Fy, Vx, Vy)
    n = prod(sz - 2); % # of interior cells
    f = zeros(n, 1); % Zero Flux Condition

    [Gx_p, L_vx, Gx_vx] = stokes1(1, sz, X, Y);
    [L_vx, Fx] = elim(1, sz, -L_vx, Fx, Vx);
    [Gx_vx, f] = elim(1, sz, Gx_vx, f,  Vx);

    [Gy_p, L_vy, Gy_vy] = stokes1(2, sz, X, Y);
    [L_vy, Fy] = elim(2, sz, -L_vy, Fy, Vy);
    [Gy_vy, f] = elim(2, sz, Gy_vy, f,  Vy);

    A = [L_vx, 0*L_vy, Gx_p; ...
         0*L_vx, L_vy, Gy_p; ...
         Gx_vx, Gy_vy, sparse(n, n)];

    f = [Fx; Fy; f];
end

function [A, f] = elim(dim, sz, A, f, V)
    dir = double((1:numel(sz)) == dim);
    I = interior(sz - dir);
    J = [find(boundary(I, -dir)); find(boundary(I, +dir))];
    f = f - A(:, J) * V(:);
    I(J) = true;

    [J1, I1] = find_boundary(I, -(~dir));
    [J2, I2] = find_boundary(I, +(~dir));
    A(:, [I1 I2]) = A(:, [I1 I2]) + A(:, [J1 J2]);
    A(:, [J; J1; J2]) = [];
end

function [J, I] = find_boundary(interior, P)
    [J, I] = boundary(interior, P);
    J = find(J);
    I = find(I);
end

function [Gp, Lv, Gv] = stokes1(dim, sz, X, Y)
    dir = double((1:numel(sz)) == dim);
    h = ones(dir + 1) / 2;
    Xs = average(X, h);
    Ys = average(Y, h);

    I = interior(sz - dir);
    Lv = laplacian(I, Xs, Ys);
    Lv = Lv(I, :);
    % full(Lv)

    K = I | shift(I, -dir);
    Gv = gradient(K, Xs, Ys, dir);
    % full(Gv)

    I = interior(sz);
    K = I & shift(I, -dir);
    Gp = gradient(K(2:end-1, 2:end-1), X(I), Y(I), dir);
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