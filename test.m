function test
clc;
e = 0.5;
A = [1 e; e 1];
f = [1 1]';
[R, g] = jacobi(A, f, 1)
v = randn(2,1);
v = apply(v, R, g, 12)
u = A \ f
