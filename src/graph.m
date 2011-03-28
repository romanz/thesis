betas = linspace(0, 10, 51);
u = zeros(size(betas));

sol = [];
for i = 1:numel(betas);
    sol = steady([], (0.1:0.1:1)*betas(i), [0 1], 2);
    u(i) = sol.Vinf;
end
