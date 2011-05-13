% SOL = SECANT(F, [A B], ITERS)
% Secant method to f(x)=0 starting at [a,b]. 
function [sol, x, y] = secant(func, x, iters, status)

    if nargin < 4
        status = @(x, y) fprintf('Secant: %e -> %e\n', x, y);
    elseif isempty(status)
        status = @(x, y) 0;
    end
    function y = f(x)
        y = func(x);
        status(x, y);
    end
    y(1) = f(x(1));
    y(2) = f(x(2));
    for k = 2:iters+1
        x(k+1) = next(x(k-1:k), y(k-1:k));
        y(k+1) = f(x(k+1));
    end
    sol = x(k+1);
end

function sol = next(x, y)
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    if dx && dy
        sol = x(1) - y(1) * dx / dy;
    else
        sol = x(2);
    end
end
