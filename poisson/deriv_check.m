syms x y
u = exp(x) + sin(y);
c = cos(x.*y);
g = [diff(u, x); diff(u, y)];
L = [diff(g(1) * c, x), diff(g(2) * c, y)];
L = simple(L);
pretty(L)
