function spheric_stokes_tb(nx, ny, iters)
x = logspace(1, 1.1, nx);
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
Ux0 = cos(gridVx.Y(:)); 
Ux0 = reshape(Ux0, gridVx.sz);
Uy0 = -sin(gridVy.Y(:)); 
Uy0 = reshape(Uy0, gridVy.sz);
P0  = zeros(gridP.sz);
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
rhsVx = rhs1x(Vx_inf) + rhs2x(Vx_inf);

[lhs1y, rhs1y] = average_dirichlet(gridVy, [-1 0]);
[lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
[lhs3y] = dirichlet(gridVy, [0 -1]);
[lhs4y] = dirichlet(gridVy, [0 +1]);
lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
rhsVy = rhs1y(Vy_inf) + rhs2y(Vy_inf);

lhs = blkdiag(lhsVx, lhsVy, speye(gridP.numel));
rhs = [rhsVx; rhsVy; sparse(gridP.numel, 1)];

% A = blkdiag(diag(gridVx.X(:).^2 .* gridVx.Y(:).^2), ...
%             diag(gridVy.X(:).^2 .* gridVy.Y(:).^2), ...
%             diag(gridP.X(:).^2 .* gridP.Y(:).^2))* A;
B = [A * lhs;];
t = [A * rhs;];

%% Timing
f = A*u;
e = f;
[Fx, Fy, div] = split(e, gridVx.sz-2, gridVy.sz-2, gridP.sz);
[Vx1, Vy1, P1] = split(u, gridVx.sz, gridVy.sz, gridP.sz);

% disp(norm(Fx(:), inf))
% disp(norm(Fy(:), inf))
% disp(norm(div(:), inf))

z0 = [Ux0(gridVx.I); Uy0(gridVy.I); P0(:)];
z = z0;
[Vx, Vy, P] = split(z, gridVx.sz-2, gridVy.sz-2, gridP.sz);
M = stokes_vanka_redblack(gridP.sz, B);
h = progress([], 0, 'Iterations');
res = NaN;
for iter = 1:iters
    for k = 1:1 % !!!
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
fprintf('Residual: %.5e\n', res);
fprintf('Vx error: %.5e\n', norm(col(Vx - Vx1(2:end-1, 2:end-1)), inf))
fprintf('Vy error: %.5e\n', norm(col(Vy - Vy1(2:end-1, 2:end-1)), inf))
fprintf('P error: %.5e\n', norm(col(P - P1), inf))

show(0, -A*u, {gridVx.sz-2, gridVy.sz-2, gridP.sz});
show(10, r, {gridVx.sz-2, gridVy.sz-2, gridP.sz});
save results

function show(f, z, sizes)
[Vx, Vy, P] = split(z, sizes{:});
mesher = @(X) mesh(X.');
figure(f+1); plot(Vx);
figure(f+2); mesher(Vy);
figure(f+3); mesher(P);
