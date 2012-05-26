function r = linreg(x, y)
    x = x(:);
    y = y(:);
    X = [x, ones(size(x))];
    h = X \ y;
    r.a = h(1);
    r.b = h(2);
    y_est = X * h;
    r.C = norm(y'*y_est) / (norm(y)*norm(y_est)); % coefficient of correlation
end