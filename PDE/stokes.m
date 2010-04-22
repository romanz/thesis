% 2D Stokes equation discretization for specific dimension
% {X,Y} are cell centers coordinates (incl. the boundary cells)
function [A, f, Mr, Mb] = stokes(sz, X, Y, Fx, Fy, Vx, Vy)
    n = prod(sz - 2); % # of interior cells
    f = zeros(n, 1); % Zero Flux Condition

    [Gx_p, L_vx, Gx_vx] = stokes1(sz, X, Y, 1);
    [L_vx, Fx] = elim(1, sz, -L_vx, Fx, Vx);
    [Gx_vx, f] = elim(1, sz, Gx_vx, f,  Vx);

    [Gy_p, L_vy, Gy_vy] = stokes1(sz, X, Y, 2);
    [L_vy, Fy] = elim(2, sz, -L_vy, Fy, Vy);
    [Gy_vy, f] = elim(2, sz, Gy_vy, f,  Vy);

    A = [L_vx, sparse(size(L_vx, 1), size(L_vy, 2)), Gx_p; ...
         sparse(size(L_vy, 1), size(L_vx, 2)), L_vy, Gy_p; ...
         Gx_vx, Gy_vy, sparse(n, n)];

    f = [Fx; Fy; f];    
    [Mr, Mb] = redblack(sz, A);    
end

function [A, f] = elim(dim, sz, A, f, V)
    dir = double((1:numel(sz)) == dim);
    I = interior(sz - dir);
    J = [find(boundary(I, -dir)); find(boundary(I, +dir))];
    f = f(:) - A(:, J) * V(:);
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

function [M1, M2] = redblack(sz, A)
    sz = sz - 2;
    K = 1:prod(sz);
    K = K(:);
    [I, J] = ind2sub(sz, K);
    
    % Same indices for variables and equations.
    V = [[index(sz - [1 0], I-1, J), index(sz - [1 0], I, J)], ...
         [index(sz - [0 1], I, J-1), index(sz - [0 1], I, J)] + prod(sz - [1 0]), ...
         K + prod(sz - [1 0]) + prod(sz - [0 1])];
    
    K = logical(mod(I - J, 2)); 
    V1 = arr2cell(V( K, :)); % odd
    V2 = arr2cell(V(~K, :)); % even
    
    M1 = vanka(A, V1, V1);
    M2 = vanka(A, V2, V2);
end

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
