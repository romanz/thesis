syms x y
u = sin(x) + cos(y)
c = exp(x - y)
g = [diff(u, x); diff(u, y)];
L = [diff(g(1) * c, x), diff(g(2) * c, y)];
L = simple(sum(L));
L
