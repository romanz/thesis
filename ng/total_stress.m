clear; load results
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
size(u)
Vx1 = [u([1, 1:n, n])'; solVx];
Vy1 = solVy;
P1 = [u(n+1:end)'; solP];

%{
P = expand([true(numel(gridVx.y), 1); false(size([Vx0(:); Vy0(:)]); zeros(numel(gridP.y), 1); P0(:)];);
q = [zeros(numel(gridVx.y), 1); Vx0(:); Vy0(:); zeros(numel(gridP.y), 1); P0(:)];

A = stokes_operator * P;
b = stokes_operator * q;
u = A \ (-b);
u = P*u + q;
[Vx, Vy, P] = split(u, gridVx.sz, gridVy.sz, gridP.sz);
%}

figure(1); mesh(Vx1(1:3, :))
figure(2); mesh(P1(1:2, :))
