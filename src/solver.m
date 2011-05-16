%% Continuation solver for a range of betas
clc; clear; tic;
betas = logspace(-2, 1, 150);
betas(2)/betas(1) - 1
sol = steady([], betas(1), [1, 0.9], 3, 'alpha', 0);
sol.alpha = 0.5;
x = 1.1;
sols = cell(numel(betas), 1);
for k = 1:numel(betas)
    sol.time = toc;
    sol = steady(sol, betas(k), sol.Vinf*[1 0.9], 3);
    sol.time = toc - sol.time;
    sols{k} = sol;
end
save
