function main3
clc
iters = 100;
x = [-1 -.3 2 3]';
y = [-2 .2 1];
Vx  = @(x, y)  y.* x.^2;
Vy  = @(x, y) -x.* y.^2;
Fx  = @(x, y) -2*y;
Fy  = @(x, y)  2*x;
div = @(x, y)  0*x + 0*y;

[sz, Xc, Yc, Xx, Yx, Xy, Yy] = create_grids(x, y);

%% Vx
[Lx, Rx, Dx, Ux, Ix, Ex] = create_boundaries(sz + [0 1]);
shave(Vx(Xx, Yx), 1, 1) % display correct Vx
Fx = Fx(Xx(Ix), Yx(Ix));
Jx = [find(Lx); find(Rx)];
Wx = ones(nnz(Jx), 1);
ux = Vx(Xx(Jx), Yx(Jx));

Dx1 = shift(Dx, [0, +1]);
Ux1 = shift(Ux, [0, -1]);
Jx = [Jx * [1 0]; ...
     [find(Dx), find(Dx1)]; ...
     [find(Ux), find(Ux1)]; ...
     find(Ex) * [1 0]]; 
Wx = [Wx * [1 0]; ...
      repmat([1, 1], nnz(Dx), 1); ...
      repmat([1, 1], nnz(Ux), 1); ...
      ones(nnz(Ex), 1) * [1 0]];
ux = [ux; ...
    Vx(Xx(Dx), Yx(Dx)) + Vx(Xx(Dx1), Yx(Dx1)); ...
    Vx(Xx(Ux), Yx(Ux)) + Vx(Xx(Ux1), Yx(Ux1)); ...
    zeros(nnz(Ex), 1)];

%% Vy
[Ly, Ry, Dy, Uy, Iy, Ey] = create_boundaries(sz + [1 0]);
shave(Vy(Xy, Yy), 1, 1) % display correct Vy
Fy = Fy(Xy(Iy), Yy(Iy));
Jy = [find(Dy); find(Uy)];
Wy = ones(nnz(Jy), 1);
uy = Vy(Xy(Jy), Yy(Jy));

Ly1 = shift(Ly, [+1, 0]);
Ry1 = shift(Ry, [-1, 0]);
Jy = [Jy * [1 0]; 
     [find(Ly), find(Ly1)]; ...
     [find(Ry), find(Ry1)]; ...
     find(Ey) * [1 0]]; 
Wy = [Wy * [1 0]; ...
      repmat([1, 1], nnz(Ly), 1); ...
      repmat([1, 1], nnz(Ry), 1); ...
      ones(nnz(Ey), 1) * [1 0]];
uy = [uy; ...
    Vy(Xy(Ly), Yy(Ly)) + Vy(Xy(Ly1), Yy(Ly1)); ...
    Vy(Xy(Ry), Yy(Ry)) + Vy(Xy(Ry1), Yy(Ry1)); ...
    zeros(nnz(Ey), 1)];

%% Build equations
[Gx_p, L_vx, Gx_vx] = stokes1(1, Xx, Yx, Xc, Yc);
[Gy_p, L_vy, Gy_vy] = stokes1(2, Xy, Yy, Xc, Yc);

[fx, L_vx ] = subst(-L_vx, Jx, Wx, ux);
[f1, Gx_vx] = subst(Gx_vx, Jx, Wx, ux);

[fy, L_vy ] = subst(-L_vy, Jy, Wy, uy);
[f2, Gy_vy] = subst(Gy_vy, Jy, Wy, uy);

% CR: sz hack!
n = prod(sz - 1); % # of interior cells
A = [L_vx, sparse(size(L_vx, 1), size(L_vy, 2)), Gx_p; ...
     sparse(size(L_vy, 1), size(L_vx, 2)), L_vy, Gy_p; ...
     Gx_vx, Gy_vy, sparse(n, n)];

f = [Fx+fx; Fy+fy; f1+f2+div(Xc(:), Yc(:))];

% CR: sz hack!
[Mr, Mb] = stokes_vanka_redblack(sz+1, A);

x = randn(numel(f), 1);
res = zeros(iters, 1);
for k = 1:iters
    r = f - A*x;    x = x + Mr*r; 
    r = f - A*x;    x = x + Mb*r;
    res(k) = norm(r, 2);
end
% CR: sz hack!
sz = sz - 1;
Vx = reshape(x(1:numel(Fx)), sz - [1 0])
Vy = reshape(x((1:numel(Fy)) + numel(Fx)), sz - [0 1])
P = reshape(x((1 + numel(Fx) + numel(Fy)):end), sz);
P = P - P(1) + 0

clf;
semilogy(res); 
axis([0 iters 1e-16 100])

save results
end

function [L, R, D, U, I, E] = create_boundaries(sz)
    I = interior(sz);
    L = ~I & shift(I, [-1, 0]);
    R = ~I & shift(I, [+1, 0]);
    D = ~I & shift(I, [0, -1]);
    U = ~I & shift(I, [0, +1]);
    E = ~(I | L | R | D | U);
end

function [sz, cX, cY, xX, xY, yX, yY] = create_grids(x, y)
    x = x(:);
    y = y(:);
    sz = [numel(x) numel(y)];

    % extended grid (for ghost points)
    x = [2*x(1) - x(2); x; 2*x(end) - x(end-1)];
    y = [2*y(1) - y(2); y; 2*y(end) - y(end-1)];

    % cell-centered coordinates
    xc = average(x, [1; 1]/2);
    yc = average(y, [1; 1]/2);

    [cX, cY] = ndgrid(xc(2:end-1), yc(2:end-1));
    [xX, xY] = ndgrid(x(2:end-1), yc);
    [yX, yY] = ndgrid(xc, y(2:end-1));
end