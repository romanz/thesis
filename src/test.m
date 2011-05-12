clc; clear;

sol = steady([], 1, [1, 2], 3, 'alpha', 0);
for a = 0:0.01:0.12
    sol = steady(sol, sol.beta, sol.Vinf*[1 1.001], 3, 'alpha', a);
    mesh(sol.C);
    disp(sol)
    drawnow;
end