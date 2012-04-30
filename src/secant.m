% Secant method solver.
% [iter, x] = secant([x0, x1])
% for i = 1:10
%   x = iter(func(x))

function [iterator, x0] = secant(x0, x1)

    n = 0;
    sol.x = [x0, x1];
    sol.y = [];
    function next_x = iter(y)
        n = n + 1;
        sol.y(n) = y;
        if n >= numel(sol.x)
            next_x = step(sol.x(n-1:n), sol.y(n-1:n));
            sol.x(n + 1) = next_x;
        else
            next_x = sol.x(n + 1);
        end
    end
    iterator = @iter;

end

function next_x = step(x, y)
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    if dx ~= 0 && dy ~= 0
        next_x = x(1) - y(1) * dx / dy;
    else
        next_x = x(2);
    end
end
