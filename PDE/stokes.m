% 2D Stokes equation discretization for specific dimension
% {X,Y} are cell centers coordinates (incl. the boundary cells)
function [A, f] = stokes(sz, X, Y, Fx, Fy, Vx, Vy)
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
