function spheric_stokes_tb(K)

nx = 30*K;
ny = 20*K;
s = 1;
x = logspace(0, .6, nx)*s;
y = linspace(0, pi, ny)*s;
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
u = [cos(gridVx.Y(:)) .* (1 - 1.5 ./ gridVx.X(:) + 0.5 ./ gridVx.X(:).^3); ...
    -sin(gridVy.Y(:)) .* (1 - 0.75 ./ gridVy.X(:) - 0.25 ./ gridVy.X(:).^3); ...
    -1.5 * cos(gridP.Y(:)) ./ gridP.X(:).^2];
% Somehow, Fx component is not O(h^2).

% u = [gridVx.X(:) .* sin(gridVx.Y(:)) .* sin(gridVx.Y(:)); ...
%      gridVy.X(:) .* cos(gridVy.Y(:)) .* sin(gridVy.Y(:)); ...
%     0*-1.5 * cos(gridP.Y(:)) ./ gridP.X(:).^2];
f = A*u;
e = f;
[Fx, Fy, div] = split(e, gridVx.sz-2, gridVy.sz-2, gridP.sz);

disp(norm(Fx(:), inf))
disp(norm(Fy(:), inf))
disp(norm(div(:), inf))
