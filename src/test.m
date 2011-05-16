%%
clc; clear;
sol = steady([], 1e-3, [1, 0.9], 3, 'alpha', 0);
sol.alpha = 0.5;
x = 1.1;
M = round(log(10e3)/log(x));
sols = cell(M, 1);
for k = 1:M
    sol.beta = (sol.beta * x);
    sol = steady(sol, sol.beta, sol.Vinf*[1 0.9], 3);
    sols{k} = sol;
end
save
%%
load
clf;
for k = 1:numel(sols)
    sol = sols{k};
    streamlines(sol, linspace(-sol.Vinf, sol.Vinf, 100));
    axis equal;
    d = 5;
    xlim([-d d]);
    ylim([0 d])
    F(k) = getframe;    
end
%%
movie(F)
