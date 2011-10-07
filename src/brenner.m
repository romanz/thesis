% Chemical Engineering Science, 1964, Vol. 19, pp. 703-727. Pergamon Press Ltd., Oxford. Printed in Great Britain.
% The Stokes resistance of an arbitrary particle-IV 
% Arbitrary fields of flow 
% H. BRENNER
% See p.711 for
function U = brenner(sol)

avg = @(x) (x(1:end-1) + x(2:end))/2;
v = [0; sol.Vs; 0];

t = sol.grid.Vy.y;
dt = diff(t);

U = -sum(avg(v .* sin(t).^2 ) .* dt) / 2;
