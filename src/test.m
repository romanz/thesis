clc; clear;

sol = steady([], 1, [1, 0.9], 3, 'alpha', 0);
for a = [0.5]
    sol = steady(sol, sol.beta, sol.Vinf*[1 0.9], 3, 'alpha', a);
    disp(sol)
end
%%
clf;
streamlines(sol, linspace(0, 1, 50));
axis equal;
xlim([-3 3]); ylim([0 3])
