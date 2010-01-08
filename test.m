function test
clc;
e = 0.5;
A = [1 e; e 1];
f = [1 1]';
[B] = jacobi(A);
v = [0 0]';

res = @(v) repmat(f, [1 size(v, 2)]) - A*v;
[v, V] = iterate(v, 0.1*B, res, 100);

u = A \ f, v
plot(sum(res(V).^2))
