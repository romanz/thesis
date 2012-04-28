% Numerical quadrature by Simpson's rule.
% >>  F = simpson(x, f);
function F = simpson(x, f)
    f = f(:);
    x = x(:);
    n = numel(x);
    assert(numel(f) == n); % same length
    assert(mod(n, 2) == 1) % odd number of points, even number of intervals
    
    [I, J] = ndgrid(1:3, 0:2:(n-2));
    K = I + J; % indexing matrix for Simpson's rule    
    w = [1;4;1]/6; % Simpson's weights, where x(i) = (x(i+1) + x(i-1))/2.
    
    f = sum(f(K) .* w(I)); % row
    dx = x(K(end, :)) - x(K(1, :)); % column
    F = f * dx;
end
