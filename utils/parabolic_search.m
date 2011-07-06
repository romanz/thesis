% Approximate f by ax+b and iterate to find the extremum of (x*f(x)).
function [result, x, y] = parabolic_search(f, x, iters)
for k = 1:2
    y(k) = f(x(k));
end
for k = 1:iters
    X = x(end-1:end);
    Y = y(end-1:end);
    dX = max(diff(X), 0.001);
    dY = diff(Y);
    a = dY / dX;
    b = (Y(1)*X(2) - Y(2)*X(1)) / dX;
    x(end+1) = -b / (2*a);
%     x(end+1) = -(Y(1)*X(2) - Y(2)*X(1)) / ((Y(2) - Y(1));
    y(end+1) = f(x(end));
    if abs(x(end) - x(end-1)) / x(end) < 1e-6 % x tolerance
        break
    end
end
result = x(end);
