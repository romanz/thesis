clc;

x0 = [1,2]
f = @(x) [x(2)-x(1); exp(-x(1)) - exp(-x(2))];
g = @(x) f(x) - f(x0);
H = @(x) [-1, 1; -exp(-x(1)), exp(-x(2))];

newton([1,0], g, H, 10)