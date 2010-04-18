function main3

sz = [6 6];
[X, Y] = ndgrid(1:sz(1), 1:sz(2));

Fx = ones(sz - [3 2]);
Fy = ones(sz - [2 3]);
[A, f, Mr, Mb] = stokes(sz, X, Y, Fx(:), Fy(:), ...
    ones(sz(2)-2, 2), ones(sz(1)-2, 2));
size(A)
x = randn(numel(f), 1);
for k = 1:1000
    x = x + Mr*(f-A*x); 
    x = x + Mb*(f-A*x);
end
Vx = reshape(x(1:numel(Fx)), sz - [3 2])
Vy = reshape(x((1:numel(Fy)) + numel(Fx)), sz - [2 3])
P = reshape(x((1 + numel(Fx) + numel(Fy)):end), sz - 2);
P = P - P(1)
