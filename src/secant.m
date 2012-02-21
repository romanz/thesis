function solver = secant(f, x)

    iter = 0;
    y = [];
    function [x0, y0] = iterator()
        iter = iter + 1;
        if iter > 2
            I = iter - (1:2);
            x(iter) = step(x(I), y(I));
        end
        x0 = x(iter);
        y0 = f(x0);
        x(iter) = x0;
        y(iter) = y0;
    end
    solver = @iterator;

end

function sol = step(x, y)
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    if dx ~= 0 && dy ~= 0
        sol = x(1) - y(1) * dx / dy;
    else
        sol = x(2);
    end
end
