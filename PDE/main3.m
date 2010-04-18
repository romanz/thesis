function main3

sz = [6 6];
[X, Y] = ndgrid(1:sz(1), 1:sz(2));

Fx = ones(sz - [3 2]);
Fy = ones(sz - [2 3]);
Vx = ones(sz(1)-2, 2);
Vy = ones(sz(2)-2, 2);
[A, f, Mr, Mb] = stokes(sz, X, Y, Fx(:), Fy(:), Vx, Vy);
size(A)
x = randn(numel(f), 1);
for k = 1:1000
    x = x + Mr*(f-A*x); 
    x = x + Mb*(f-A*x);
end
% x(1:(numel(Fx)+numel(Fy)))
x