clc
syms x y
u = x*y
c = x*y;
g = [diff(u, x); diff(u, y)]
L = [diff(g(1) * c, x), diff(g(2) * c, y)];
L = simple(sum(L))
v = [x+y, x-y];
flux = diff(v(1), x) + diff(v(2), y)
A = simple(v * g)
S = simple(L - A)
