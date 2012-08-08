clc; 
clear;
D = dir('for_plot/*.mat');

for k = 1:numel(D)
    d = D(k);
    f = d.name;
    S{k} = load(['for_plot/', f]);
end

S = reshape(S, 3, numel(S)/3);
S = S(:, 2:end);

betas = cellfun(@(s) s.betas(end), S);
[betas, I] = sort(betas, 2);
V = cellfun(@(s) s.V(end), S);

Vn = V(:, I(1, :));
Vr0 = [-1/3, 4/3] * Vn(1:2, :);
Vr1 = [-1/3, 4/3] * Vn(2:3, :);
Vr2 = Vr1 + (Vr1 - Vr0)/15;

U1 = 4.25752957407993;
U3 = 0.966869;

b = betas(1, :);

V1 = U1*b;
V3 = U3*b.^3;

figure(1); loglog(b, Vr1-V1, 'o-', b, V3, 's-')
legend('Numerical - Linear', 'Analytical cubic correction', ...
    'Location', 'SouthEast')
grid on
xlabel('\beta');
ylabel('steady-state velocity: deviation from linear')
title('Numerical results for cubic correction')
print -depsc2 graph

%{
figure(2); loglog(b, Vr1, 'o-', b, V1+V3, 's-', b, V1, 'd-')
legend('Numerical', 'Cubic', 'Linear', 'Location', 'SouthEast')
grid on
xlabel('\beta'); ylabel('steady-state velocity')
%}