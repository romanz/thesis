function S = total_stress(solVx, solVy, solP, radius, theta)
Vx0 = solVx(1:2, :);
Vy0 = solVy(1:2, :);
P0 = solP(1, :);
[~, gridP, gridVx, gridVy] = ...
    grids([2 -1; 1 0; 0 1]*radius(1:2), theta);
stokes_operator = spheric_stokes(gridVx, gridVy, gridP);        
n = numel(gridP.y);
% Note correct axis choice!
I = [true(nnz(gridVx.I), 1); false(nnz(gridVy.I), 1); repmat([true; false], n, 1)];
stokes_operator1 = stokes_operator(I, :);
% Note correct axis choice!
J = [col([0 0 0; repmat([1 0 0], n, 1); 0 0 0]'); false(gridVy.numel, 1); repmat([1; 0], n, 1)];
A = stokes_operator1 * expand(J);
% Note correct axis choice!
q = stokes_operator1 * [col([zeros(1, size(Vx0, 2)); Vx0]); col([zeros(1, size(Vy0, 2)); Vy0; zeros(1, size(Vy0, 2))]); col([0*P0(:) P0(:)]')];
u = A \ -q;
Vx0 = [u([1, 1:n, n])'; Vx0];
P0 = [u(n+1:end)'; P0];


Sr = -mean(P0) ...
     +diff(Vx0([1 3], 2:end-1)) / diff(radius(1:2));
St =  average(diff(Vy0) / diff(radius(1:2)), [1 1]/2) ...
     -average(Vy0(1:2, :), [1 1; 1 1]/4);
t = gridP.y';
F = -Sr .* cos(t) + St .* sin(t);
S = sum((2*pi*sin(gridP.y)) .* F' .* diff(theta));    
