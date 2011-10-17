clear;
syms a r t U W
clc
at = @(f) simple(subs(f, r, 1));

% psi = -(U/2) * (r^2 - (3*a*r)/2 + a^3/(2*r))*sin(t)^2;
% V = curl(psi);
% L = simple(vector_laplacian(V));
% P = -3*U*a*cos(t)/(2*r^2);
% gP = gradient(P);
% assert(all(L == gP))
% F = radial_stress(V, P, a);

psi = -(W/2) * (r^2 - a^3/r)*sin(t)^2;
V = curl(psi);
assert([1 0]*subs(V, r, a) == 0);
L = simple(vector_laplacian(V));
L1 = vlapl(V);
eL = (L - L1)

Phi1 = ((1/4)*r^-2 - r) * cos(t);
C1 = (3/4)*r^-2 * cos(t);
assert(subs(Phi1 + C1, r, 1) == 0)
assert(subs(diff(Phi1 - C1, r), r, 1) == 0)
syms gamma phi
xi = -phi-log(gamma);
g = 4*log((exp(xi/2)+1)/2);
V1a = subs(g, phi, 0) * at(diff(Phi1, t));
U1 = simple(int(V1a * sin(t)^2, 0, pi) / int(sin(t), 0, pi));

V1 = subs(V, a, 1);
W1 = simple(int(at(V1(2)) * sin(t)^2, 0, pi) / int(sin(t), 0, pi));
assert(W1 == W)

gC1 = gradient(C1);
gPhi1 = gradient(Phi1);

syms A1 A2 A3 A4

RHS_Phi = simple(-gC1.'*gPhi1);
Phi2 = 3/32*r^-4*sin(t)^2 - 3/8*r^-1*sin(t)^2 - 3/32*r^-4;
Phi2 = simple(Phi2);
assert(simple(scalar_laplacian(Phi2) - RHS_Phi) == 0)

RHS_C = simple(V1.'*gC1);
C2 = 3*W/8*((r^-1 + 1/2*r^-4)*sin(t)^2 - 1/2*r^-4);
C2 = simple(C2);
s = simple(RHS_C - scalar_laplacian(C2));
assert(s == 0);

C2 = C2 + A1/r + A2*(3*cos(t)^2-1)/r^3;
Phi2 = Phi2 + A3/r + A4*(3*cos(t)^2-1)/r^3;

Dr = @(f) diff(f, r);
s1 = at(C2 + Phi2 - Phi1^2/2);
s2 = at(Dr(C2 - Phi2 - C1*Dr(Phi1)));

syms pi
q1 = subs(s1, t, 0);
q2 = subs(s1, t, pi/2);
q3 = subs(s2, t, 0);
q4 = subs(s2, t, pi/2);

S = solve(q1, q2, q3, q4, A1, A2, A3, A4);
C2 = subs(C2, [A1 A2 A3 A4], [S.A1 S.A2 S.A3 S.A4]);
Phi2 = subs(Phi2, [A1 A2 A3 A4], [S.A1 S.A2 S.A3 S.A4]);

assert(simple(scalar_laplacian(Phi2) - RHS_Phi) == 0)
assert(simple(RHS_C - scalar_laplacian(C2)) == 0)

assert( 0 == at(C2 + Phi2 - Phi1^2/2) );
assert( 0 == at(Dr(C2 - Phi2 - C1*Dr(Phi1))) );

V2 = subs(g, phi, 0) * at(diff(Phi2, t)) + ...
    subs(diff(g, phi), phi, 0) * at(Phi1 * diff(Phi1, t));
V2 = simple(V2);

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
