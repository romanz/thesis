% function diffusion_advection
clear
sz = [2^3+1, 1];
V = -2;
N = prod(sz);
x = linspace(0, 1, sz(1));
y = 0;
[X, Y] = ndgrid(x, y);
C = 1 + 0*X + 0*Y;
[L, M, I] = laplacian(sz, X, Y, C);
L = dinv(M) * L;
[Gx, Gy] = gradient(sz, X, Y);
A = L - V*Gx - 0*Gy;
f = zeros(sz);

Bl = boundary(I, [-1 0]); % Left boundary (small X)
Br = boundary(I, [+1 0]); % Right boundary (large X)
[A, f] = dirichlet(A, f, Bl, X, Y, @(X,Y) 0*X+0*Y+1);
[A, f] = dirichlet(A, f, Br, X, Y, @(X,Y) 0*X+0*Y);
A = A(I, I); f = f(I);
t = linspace(min(x), max(x), 10e3);
u = A\f;
plot(x, [1; u; 0], t, 1-(exp(V*t)-1)/(exp(V)-1))

