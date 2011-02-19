% -Laplacian{Vx,Vy} + Grad{P} = Operator{Vx,Vy,P}
function A = spheric_stokes(gridVx, gridVy, gridP)

[Gp_x, Lv_x0, Dv_x] = stokes1(1, gridVx, gridP);
[Gp_y, Lv_y0, Dv_y] = stokes1(2, gridVy, gridP);

Ix = find(gridVx.I);
Iy = find(gridVy.I);

% -2Vr / r^2
Lv_x1 = sparse(1:numel(Ix), Ix, ...
    -2./gridVx.X(Ix).^2, numel(Ix), gridVx.numel);

% -2/(r^2 sin(theta)) d(Vtheta sin(theta))/dtheta
Lv_x2 = sparse(1:numel(Ix), 1:numel(Ix), ...
    -2./( sin(gridVx.Y(Ix)) .* gridVx.X(Ix).^2 )) * ...
    expand([false(1, gridVx.sz(2)-2); true(gridVx.sz - 2); false(1, gridVx.sz(2)-2)])' * ...
    interpolator(average(gridVy.X, [1, 1]/2), gridVx.X(:, 2:end-1)) * ...
    deriv(gridVy.Y, [0 1]) * ...
    sparse(1:prod(gridVy.sz), 1:prod(gridVy.sz), sin(gridVy.Y));

% -Vtheta/(r * sin theta)^2
Lv_y1 = -sparse(1:numel(Iy), Iy, ...
    1./( gridVy.X(Iy) .* sin(gridVy.Y(Iy)) ).^2, numel(Iy), gridVy.numel);

% 2/r^2 dVr/dtheta
Lv_y2 = sparse(1:numel(Iy), 1:numel(Iy), 2./gridVy.X(Iy).^2) * ...
    expand([false(gridVy.sz(1)-2, 1), true(gridVy.sz - 2), false(gridVy.sz(1)-2, 1)])' * ...
    interpolator(average(gridVx.X, [1, 1]/2), gridVy.X(2:end-1, :)) * ...
    deriv(gridVx.Y, [0 1]);

Z = sparse(prod(gridP.sz), prod(gridP.sz));

A = [-Lv_x0 - Lv_x1, -Lv_x2,         Gp_x; ...
     -Lv_y2,         -Lv_y0 - Lv_y1, Gp_y; ...
     Dv_x,            Dv_y,          Z];

