nx = 160;
ny = 100;
s = 1;
x = logspace(0, 1, nx)*s;
y = linspace(0, pi, ny)*s;
[center, gridP, gridVx, gridVy] = grids(x, y);
A = spheric_stokes(gridVx, gridVy, gridP);
u = [1*gridVx.X(:); 0*gridVy.X(:); 0*gridP.X(:)];
f = A*u;
[Fx, Fy, div] = split(f, gridVx.sz-2, gridVy.sz-2, gridP.sz);
Fx
norm(Fx(:), inf)
norm(Fy(:), inf)
norm(div(:), inf)
