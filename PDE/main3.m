function main3
clc
iters = 40;
sz = [1 1]*2^2;
[X, Y] = ndgrid(1:sz(1), 1:sz(2));

% Vx
Vx = @(X, Y) 1 + 0*X.^2 .* Y;
[Xx, Yx, Lx, Rx, Dx, Ux, Ix, Ex] = staggered_boundary(sz, X, Y, 1);
Fx = @(x, y) 0*x+0*y+0;
Fx = Fx(Xx(Ix), Yx(Ix));
Jx = [find(Lx); find(Rx)];
Wx = 1+0*Jx;
ux = Vx(Xx(Jx), Yx(Jx));

Jx = [Jx * [1 0]; ...
     [find(Dx), find(shift(Dx, [0, +1]))]; ...
     [find(Ux), find(shift(Ux, [0, -1]))]; find(Ex) * [1 0]]; 
Wx = [Wx * [1 0]; repmat([-1, 1], nnz(Dx), 1); repmat([1, -1], nnz(Ux), 1); ...
      ones(nnz(Ex), 2)];
ux = [ux; 0*find(Dx); 0*find(Ux); zeros(nnz(Ex), 1)];

% Vy
Vy = @(X, Y) 1 - 0*Y.^2 .* X;
[Xy, Yy, Ly, Ry, Dy, Uy, Iy, Ey] = staggered_boundary(sz, X, Y, 2);
Fy = @(x, y) 0*x+0*y+0;
Fy = Fy(Xy(Iy), Yy(Iy));
Jy = [find(Dy); find(Uy)];
Wy = 1+0*Jy;
uy = Vy(Xy(Jy), Yy(Jy));

Jy = [Jy * [1 0]; ...
     [find(Ly), find(shift(Ly, [+1, 0]))]; ...
     [find(Ry), find(shift(Ry, [-1, 0]))]; find(Ey) * [1 0]]; 
Wy = [Wy * [1 0]; repmat([-1, 1], nnz(Ly), 1); repmat([1, -1], nnz(Ry), 1); ...
      ones(nnz(Ey), 2)];
uy = [uy; 0*find(Ly); 0*find(Ry); zeros(nnz(Ey), 1)];

[Gx_p, L_vx, Gx_vx] = stokes1(sz, X, Y, 1);
[Gy_p, L_vy, Gy_vy] = stokes1(sz, X, Y, 2);

n = prod(sz - 2); % # of interior cells

[fx, L_vx ] = subst(-L_vx, Jx, Wx, ux);
[f1, Gx_vx] = subst(Gx_vx, Jx, Wx, ux);

[fy, L_vy ] = subst(-L_vy, Jy, Wy, uy);
[f2, Gy_vy] = subst(Gy_vy, Jy, Wy, uy);

A = [L_vx, sparse(size(L_vx, 1), size(L_vy, 2)), Gx_p; ...
     sparse(size(L_vy, 1), size(L_vx, 2)), L_vy, Gy_p; ...
     Gx_vx, Gy_vy, sparse(n, n)];

f = [Fx+fx; Fy+fy; f1+f2];    

%
[Mr, Mb] = stokes_vanka_redblack(sz, A);

x = randn(numel(f), 1);
res = zeros(iters, 1);
for k = 1:iters
    r = f - A*x;
    x = x + Mr*r; 

    r = f - A*x;
    x = x + Mb*r;

    res(k) = norm(r, 2);
end
Vx = reshape(x(1:numel(Fx)), sz - [3 2])
Vy = reshape(x((1:numel(Fy)) + numel(Fx)), sz - [2 3])
P = reshape(x((1 + numel(Fx) + numel(Fy)):end), sz - 2);
P = P - P(1)

clf;
semilogy(res); 
axis([0 iters 1e-16 100])

save results
end

function [X, Y, L, R, D, U, I, E] = staggered_boundary(sz, X, Y, dim)
    [X, Y, P] = staggered_grid(sz, X, Y, dim);
    I = interior(sz - P);
    L = ~I & shift(I, [-1, 0]);
    R = ~I & shift(I, [+1, 0]);
    D = ~I & shift(I, [0, -1]);
    U = ~I & shift(I, [0, +1]);
    E = ~(I | L | R | D | U);
end
