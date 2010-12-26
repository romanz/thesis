% Newton method solver.
function [x, history] = newton(x, func, hessian, iters)
    history = zeros(numel(x), iters);
    for i = 1:iters
        f = func(x);
        H = hessian(x);
        dx = -H \ f;
        history(:, i) = x;
        x = x + reshape(dx, size(x));
    end
end