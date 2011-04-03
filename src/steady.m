% Steady-state solver
% Given Beta, find Velocity for zero particle force.
function sol = steady(sol, betas, v, iters, varargin)

conf = {'version', 0, 'logging', 0, varargin{:}};
function f = func(u)
    sol = main(sol, betas(end), u, conf{:});
    f = sol.force.total;
end

v = [min(v), max(v)];
sol = main(sol, betas, v(1), conf{:});
status = @(x, f) fprintf('B = %e, V = %e, F = %e\n', betas(end), x, f);
secant(@func, v, iters, status);

end
