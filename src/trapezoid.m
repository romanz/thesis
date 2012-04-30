% Numerical quadrature by Trapezoid's rule.
% >>  F = simpson(x, f);
function F = trapezoid(x, f)
    f = f(:);
    x = x(:);
    n = numel(x);
    assert(numel(f) == n); % same length
    f = (f(1:end-1) + f(2:end)) / 2;
    dx = diff(x);
    F = f' * dx;
end
