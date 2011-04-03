betas = logspace(-1, 0.7, 50);
u = zeros(size(betas));
sol = steady([], 0, [0.1, 0.0], 2);
for i = 1:numel(betas);
    sol = steady(sol, betas(i), sol.Vinf + [0 1], 2);
    u(i) = sol.Vinf;
    fprintf('\tb = %e\n\tu = %e\n', betas(i), u(i));
end
