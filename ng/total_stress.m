function S = total_stress(solVx, solVy, solP, radius, theta)

Vx0 = solVx(1:2, :); % radial components
Vy0 = solVy(1:2, :); % tangential components
P0 = solP(1, :); 
% grid for computation of (V,P) around R=1.
[~, gridP, gridVx, gridVy] = grids([2 -1; 1 0; 0 1]*radius(1:2), theta);
stokes_operator = spheric_stokes(gridVx, gridVy, gridP);
n = numel(gridP.y); % number of cells on R=1.

% Equations extraction: 
% * include radial forces and divergence for r<1.
% * exclude tangential forces and divergence for r>1.
I = [true(nnz(gridVx.I), 1); false(nnz(gridVy.I), 1); repmat([true; false], n, 1)];
stokes_operator1 = stokes_operator(I, :);

% Variables extraction:
% * The unknowns are radial velocities and pressure for r<1.
J = [col([0 0 0; repmat([1 0 0], n, 1); 0 0 0]'); false(gridVy.numel, 1); repmat([1; 0], n, 1)];
A = stokes_operator1 * expand(J);
% * All other variables have been computed: 
%   radial velocities (r>=1), pressures (r>1)
q = stokes_operator1 * [col([zeros(1, size(Vx0, 2)); Vx0]); ...
                        col([zeros(1, size(Vy0, 2)); Vy0; zeros(1, size(Vy0, 2))]); ...
                        col([0*P0(:) P0(:)]')];
% Solve A*u + q = 0:
u = A \ -q;
% Expand radial velocity and pressure for r<1.
Vx0 = [u([1, 1:n, n])'; Vx0];
P0 = [u(n+1:end)'; P0];

% Radial force: -P + 2 dVr/dr
Fr = -mean(P0) ...
     +2*diff(Vx0([1 3], 2:end-1)) / diff(gridVx.x([1 3]));
% Tangential force: dVt/dr + 1/r {dVr/dt - Vt}
% NOTE: Vr=0 for r=1.
Ft =  average(diff(Vy0) / diff(radius(1:2)), [1 1]/2) ...
     -average(Vy0(1:2, :), [1 1; 1 1]/4);

% Quadrature (midpoint)
a = gridP.y'; % Angle axis (from 0 to pi).
F = -Fr .* cos(a) + Ft .* sin(a); % Projection on axis of symmetry.
S = sum( F' .* (2*pi*sin(gridP.y)) .* diff(theta));    

end