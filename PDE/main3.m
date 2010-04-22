function main3
clc
sz = [1 1]*4;
[X, Y] = ndgrid(1:sz(1), 1:sz(2));

% Vx
[Xx, Yx, Lx, Rx, Dx, Ux, Ix] = staggered_boundary(sz, X, Y, 1);
Fx = @(x, y) 0*x+0*y+1;
Fx = Fx(Xx(Ix), Yx(Ix));

% Vy
[Xy, Yy, Ly, Ry, Dy, Uy, Iy] = staggered_boundary(sz, X, Y, 2);
Fy = @(x, y) 0*x+0*y+1;
Fy = Fy(Xy(Iy), Yy(Iy));

[A, f] = stokes(sz, X, Y, Fx, Fy, ...
    ones(sz(2)-2, 2), ones(sz(1)-2, 2));
[Mr, Mb] = stokes_vanka_redblack(sz, A);

% size(A)
x = randn(numel(f), 1);
for k = 1:100
    x = x + Mr*(f-A*x); 
    x = x + Mb*(f-A*x);
end
Vx = reshape(x(1:numel(Fx)), sz - [3 2])
Vy = reshape(x((1:numel(Fy)) + numel(Fx)), sz - [2 3])
P = reshape(x((1 + numel(Fx) + numel(Fy)):end), sz - 2);
P = P - P(1)

norm(f - A*x, inf)

save results
end

function [X, Y, L, R, D, U, I] = staggered_boundary(sz, X, Y, dim)
    [X, Y, P] = staggered_grid(sz, X, Y, dim);
    I = interior(sz - P);
    L = ~I & shift(I, [-1, 0]);
    R = ~I & shift(I, [+1, 0]);
    D = ~I & shift(I, [0, -1]);
    U = ~I & shift(I, [0, +1]);
end
