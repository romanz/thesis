% function diffusion_advection
clear
sz = [2^5+1, 1];
V = 32;
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

Bl = ~I & circshift(I, [-1 0]); % Left boundary (small X)
Br = ~I & circshift(I, [+1 0]); % Right boundary (large X)
[A, f] = boundary_dirichlet(A, f, Bl, X, Y, @(X,Y) 0*X+0*Y);
[A, f] = boundary_dirichlet(A, f, Br, X, Y, @(X,Y) 0*X+0*Y+1);
A = A(I, I); f = f(I);
t = x(I);
plot(t, A\f, t, (exp(V*t)-1)/(exp(V)-1))
