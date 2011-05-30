% Sherman–Morrison–Woodbury formula
function test_woodbury
n = 5;
A = randn(n);
u = randn(n, 1);
v = randn(1, n);

invA = inv(A);
% (A + uv)*x = I
% [A  u][x] = [I]
% [v -1][y] = [0]
% Ax + uy = I
% x + invA uy = invA
% vx + v invA uy = v invA
% y (1 + v invA u) = v invA
% x = invA (I + u y)
% x = invA + invA u v invA / (1 + v invA u)
y = v*invA / (1 + v*invA*u);
x = invA * (eye(n) - u*y);
x - inv(A + u*v)
