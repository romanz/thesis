function show1
close all;
U1 = 4.25753;

S{1} = get('[100x19]_1e+07_alpha=0.00.mat');
S{2} = get('[200x39]_1e+07_alpha=0.00.mat');
S{3} = get('[400x79]_1e+07_alpha=0.00.mat');
U3 = 0.0178949;
plot1(1, S, U1, U3)

S{1} = get('[100x19]_1e+07_alpha=0.25.mat');
S{2} = get('[200x39]_1e+07_alpha=0.25.mat');
S{3} = get('[400x79]_1e+07_alpha=0.25.mat');
U3 = 0.357798;
plot1(2, S, U1, U3)

S{1} = get('[100x19]_1e+07_alpha=0.50.mat');
S{2} = get('[200x39]_1e+07_alpha=0.50.mat');
S{3} = get('[400x79]_1e+07_alpha=0.50.mat');
U3 = 0.966869;
plot1(3, S, U1, U3)

function [s] = get(f)
s = load(f, 'betas', 'V', 'alpha');

function plot1(fig, S, U1, U3)
figure(fig);
b = S{1}.betas;
a = S{1}.alpha;
plot(b, S{1}.V, '.', b, S{2}.V, '.', b, S{3}.V, '.', ...
    b, U1*b, b, U1*b + U3*b.^3)
xlabel('\beta')
ylabel('Steady-state velocity')
title(sprintf('\\alpha = %.2f', a))

legend('[100x20]', '[200x40]', '[400x80]', 'Linear', 'Linear + Cubic', ...
    'Location', 'SouthEast')

print('-depsc', sprintf('%d.eps', fig))