syms a r t U W
clc
psi = (U/2) * (r^2 - (3*a*r)/2 + a^3/(2*r))*sin(t)^2;
V = curl(psi);
L = simple(vector_laplacian(V));
P = -3*U*a*cos(t)/(2*r^2);
gP = gradient(P);
assert(all(L == gP))
F = radial_stress(V, P, a);

psi = (W/2) * (a*r - a^3/r)*sin(t)^2;
V = curl(psi);
assert([1 0]*subs(V, r, a) == 0);
L = simple(vector_laplacian(V));
L1 = vlapl(V);
eL = (L - L1)

Phi1 = ((1/4)*r^-2 - r) * cos(t);
C1 = (3/4)*r^-2 * cos(t);
assert(subs(Phi1 + C1, r, 1) == 0)
assert(subs(diff(Phi1 - C1, r), r, 1) == 0)
V1 = subs(V, a, 1);

gC1 = gradient(C1);
gPhi1 = gradient(Phi1);

RHS1 = simple(-gC1.'*gPhi1);
RHS2 = simple(V1.'*gC1);
return
P = W*a*cos(t)/r^2;
gP = gradient(P);
assert(all(L == gP))
F = radial_stress(V, P, a);

psi1 = (U/2) * (r^2 - (3*a*r)/2 + a^3/(2*r))*sin(t)^2;
psi2 = (W/2) * (a*r - a^3/r)*sin(t)^2;
psi2 = subs(psi2, W, 3*U/2);
psi = simple(psi1 + psi2);
pretty(psi)
V = simple(curl(psi))
pretty(V)
L = simple(vector_laplacian(V))    
P = 0; % Note that the pressure is constant!
gP = gradient(P);
assert(all(L == gP))
F = radial_stress(V, P, a)
save results
