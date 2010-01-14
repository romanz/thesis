function test
clc;
N = 2^6+1;
randn('state', 0);
u = randn(N, 1);
A = matr(N);
f = A*u;
v = zeros(N, 1);

for i = 1:50
    v = V_cycle(@matr, v, f, struct('n1', 2, 'n2', 1));
end
% u = A \ f; v;
% semilogy(sum(residual(V).^2))
plot([matr(N)\f v])

function [A, B] = matr(N)
A = tridi(-1, repmat(2, [N 1]), -1)/4;
B = jacobi(A)/3;

function v = V_cycle(matr, v, f, conf)
N = numel(v);
[A, B] = matr(N);
if N == 2
    v = A \ f;
    return;
end
residual = @(v) f - A*v;
v = iterate(v, B, residual, conf.n1);
r = residual(v);
r = r(1:2:end);
e = V_cycle(matr, zeros((N+1)/2, 1), r, conf);
e = interp1(1:2:N, e, 1:N);
v = v + e(:);
v = iterate(v, B, residual, conf.n2);
