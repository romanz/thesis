function [a, b, e, R] = linreg(x, y)
X = [x ones(size(x, 1), 1)];
w = X \ y;
a = w(1:end-1);
b = w(end);
e = y - X*w;
R = norm(e) / norm(y);
