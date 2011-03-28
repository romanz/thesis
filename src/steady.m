function sol = steady(sol, betas, v, iters)

conf = {'version', 0};
function f = func(u)
    sol = main(sol, betas(end), u, conf{:});
    f = sol.force.total;
end

v = [min(v), max(v)];
sol = main(sol, betas, mean(v), conf{:});
secant(@func, v, iters);

end
