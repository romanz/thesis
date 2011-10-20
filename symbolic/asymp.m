clear;
syms b r t A
c1 = 1 + b * 3/4 * r^(-2) * cos(t) + b^2*A;
phi1 = b * (1/4 * r^(-2) - r) * cos(t);
Dr = @(f) diff(f, r);
Dt = @(f) diff(f, t);

cond = @(phi, c) subs([phi + log(c); c*Dr(phi) - Dr(c)], r, 1);
taylor( cond(phi1, c1), 3, b )
