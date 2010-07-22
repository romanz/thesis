% >> spheric_stokes_tb(400, 30, 100e3);

function spheric_stokes_tb(nx, ny, iters)
x = logspace(0, 2, nx);
y = linspace(0, pi, ny);
% y = ( 1 - cos(linspace(0, pi, ny)) )*pi/2;
[center, gridP, gridVx, gridVy] = grids(x, y);
A = spheric_stokes(gridVx, gridVy, gridP);

% %%
% u = [1*gridVx.X(:); 0*gridVy.X(:); 0*gridP.X(:)];
% f = A*u;
% [Fx, Fy, div] = split(f, gridVx.sz-2, gridVy.sz-2, gridP.sz);
% 
% norm(Fx(:), inf)
% norm(Fy(:), inf)
% norm(div(:), inf)

% %%
% u = [gridVx.X(:).^2; 0*gridVy.X(:); 0*gridP.X(:)];
% f = A*u;
% [Fx, Fy, div] = split(f, gridVx.sz-2, gridVy.sz-2, gridP.sz);
% 
% norm(Fx(:), inf)
% norm(Fy(:), inf)
% norm(div(:), inf)

%%
Ux0 = cos(gridVx.Y(:)) .* (1 - 1.5 ./ gridVx.X(:) + 0.5 ./ gridVx.X(:).^3); 
Ux0 = reshape(Ux0, gridVx.sz);
Uy0 = -sin(gridVy.Y(:)) .* (1 - 0.75 ./ gridVy.X(:) - 0.25 ./ gridVy.X(:).^3); 
Uy0 = reshape(Uy0, gridVy.sz);
P0  = -1.5 * cos(gridP.Y(:)) ./ gridP.X(:).^2;
u = [Ux0(:); Uy0(:); P0(:)];

% u = [gridVx.X(:) .* sin(gridVx.Y(:)) .* sin(gridVx.Y(:)); ...
%      gridVy.X(:) .* cos(gridVy.Y(:)) .* sin(gridVy.Y(:)); ...
%     0*-1.5 * cos(gridP.Y(:)) ./ gridP.X(:).^2];
% size(A)
vel = 1;
Vx_inf = Ux0(end, 2:end-1); vel*cos(gridVx.y(2:end-1));
Vy_inf = Uy0(end, 2:end-1); -vel*sin(gridVy.y(2:end-1));

[lhs1x, rhs1x] = dirichlet(gridVx, [-1 0]);
[lhs2x, rhs2x] = dirichlet(gridVx, [+1 0]);
[lhs3x] = neumann(gridVx, [0 -1]);
[lhs4x] = neumann(gridVx, [0 +1]);
lhsVx = lhs1x * lhs2x * lhs3x * lhs4x * restrict(gridVx.I);
rhsVx = rhs1x(0) + rhs2x(Vx_inf);

[lhs1y, rhs1y] = average_dirichlet(gridVy, [-1 0]);
[lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
[lhs3y] = dirichlet(gridVy, [0 -1]);
[lhs4y] = dirichlet(gridVy, [0 +1]);
lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
rhsVy = rhs1y(0) + rhs2y(Vy_inf);

lhs = blkdiag(lhsVx, lhsVy, speye(gridP.numel));
rhs = [rhsVx; rhsVy; sparse(gridP.numel, 1)];

A = [A];
B = [A * lhs;];
t = [A * rhs;];

%% Timing
tic;
f = A*u;
toc;
e = f;
[Fx, Fy, div] = split(e, gridVx.sz-2, gridVy.sz-2, gridP.sz);
[Vx1, Vy1, P1] = split(u, gridVx.sz, gridVy.sz, gridP.sz);

% disp(norm(Fx(:), inf))
% disp(norm(Fy(:), inf))
% disp(norm(div(:), inf))

z = zeros(numel(t), 1);
M = stokes_vanka_redblack(gridP.sz, B);
h = progress([], 0, 'Iterative Solver');
for iter = 1:iters
    for k = 1:2
        r = t - B*z;
        z = z + M{k}*r;        
    end
    [Vx, Vy, P] = ...
        split(z, gridVx.sz-2, gridVy.sz-2, gridP.sz);
    z = [Vx(:); Vy(:); P(:)-mean(P(:))];
    res = norm(r, inf);
    progress(h, iter/iters, ...
        sprintf('Residual = %.3e', res));
end
progress(h, []);

norm(col(P - P1), inf)
norm(col(Vx - Vx1(2:end-1, 2:end-1)), inf)
norm(col(Vy - Vy1(2:end-1, 2:end-1)), inf)

save results
