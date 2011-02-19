g = init_grid(logspace(0, 5, 200), linspace(0, pi, 101));
C = 10*rand(g.sz);
L = laplacian(g.I, g.X, g.Y, C);
[D1, G1, I1] = operators(g, 1);
[D2, G2, I2] = operators(g, 2);
Lnew = D1 * spdiag(I1 * C(:)) * G1 + D2 * spdiag(I2 * C(:)) * G2;
E = abs(L - Lnew);
max(E(:))