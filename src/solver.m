%% Continuation solver for a range of betas
clc; clear;
betas = logspace(0, 1, 20);
sol = steady([], betas(1), [1, 0.9], 3, 'alpha', 0);
sol.alpha = 0.5;
sols = cell(numel(betas), 1);
tic;
for k = 1:numel(betas)
    sol.time = toc;
    sol = steady(sol, betas(k), sol.Vinf*[1 0.9], 3);
    sol.time = toc - sol.time;
    sols{k} = sol;
end
toc
save
