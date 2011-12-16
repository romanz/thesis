syms e
g = 1 + e;
U1 = 2*log((1 + g.^(-0.5))/2);
U3 = 31./(320*(sqrt(g) + 1)) - 9./(320*(sqrt(g) + 1).^2) + 1/1680 - U1*11/320;
b = sqrt(-U1./U3);

e1 = logspace(-5, 5, 100);
b1 = subs(b, e, e1);
loglog(e1, b1.^2, e1, (13440/1129)*e1, e1, 5.36^2+0*e1); 
grid