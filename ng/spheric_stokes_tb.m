function spheric_stokes_tb(nx, ny, iters)
x = logspace(0, 1, nx);
y = linspace(0, pi, ny);
[center, gridP, gridVx, gridVy] = grids(x, y);
A = spheric_stokes(gridVx, gridVy, gridP);


%%
Ux0 = cos(gridVx.Y) .* (1 - 1.5 ./ gridVx.X + 0.5 ./ gridVx.X.^3); 
Uy0 = -sin(gridVy.Y) .* (1 - 0.75 ./ gridVy.X - 0.25 ./ gridVy.X.^3); 
P0  = -1.5 * cos(gridP.Y) ./ gridP.X.^2;
u = [Ux0(:); Uy0(:); P0(:)];

[lhs1x, rhs1x] = dirichlet(gridVx, [-1 0]);
[lhs2x, rhs2x] = dirichlet(gridVx, [+1 0]);
[lhs3x] = neumann(gridVx, [0 -1]);
[lhs4x] = neumann(gridVx, [0 +1]);
lhsVx = lhs1x * lhs2x * lhs3x * lhs4x * restrict(gridVx.I);
rhsVx = rhs1x(0) + rhs2x(Ux0(end, 2:end-1));

[lhs1y, rhs1y] = average_dirichlet(gridVy, [-1 0]);
[lhs2y, rhs2y] = dirichlet(gridVy, [+1 0]);
[lhs3y] = dirichlet(gridVy, [0 -1]);
[lhs4y] = dirichlet(gridVy, [0 +1]);
lhsVy = lhs1y * lhs2y * lhs3y * lhs4y * restrict(gridVy.I);
rhsVy = rhs1y(0) + rhs2y(Uy0(end, 2:end-1));

lhs = blkdiag(lhsVx, lhsVy, speye(gridP.numel));
rhs = [rhsVx; rhsVy; sparse(gridP.numel, 1)];

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
z = randn(size(z0));
[Vx, Vy, P] = split(z, gridVx.sz-2, gridVy.sz-2, gridP.sz);
M = stokes_vanka_redblack(gridP.sz, B);
h = progress([], 0, 'Iterations');
r = t - B*z;
res = norm(r, inf);
for iter = 1:iters
    for k = 1:numel(M) 
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

show(0,  z, {gridVx.sz-2, gridVy.sz-2, gridP.sz});
show(10, z-z0, {gridVx.sz-2, gridVy.sz-2, gridP.sz});
save results

function show(f, z, sizes)
[Vx, Vy, P] = split(z, sizes{:});
mesher = @(x, y, Z) mesh(x, 180*y/pi, Z.');
figure(f+1); mesher(gridVx.x(2:end-1), gridVx.y(2:end-1), Vx);
figure(f+2); mesher(gridVy.x(2:end-1), gridVy.y(2:end-1), Vy);
figure(f+3); mesher(gridP.x, gridP.y, P);
end

end