function main_show(matname)

clf
fprintf('[%s]', matname)
load(matname)

m = 1;
I = 1:m:numel(betas);
betas = betas(I);
V = V(I);

% Linear regime steady-state velocity.
U = betas*(sol.Du*log(16) + sol.zeta)/(1 + 2*sol.Du); 

% Show numerical results and the error from linear regime
subplot 121
loglog(betas, V, '.', betas, abs(U-V), '-'); 
t = sprintf(['Electrophoresis\nDu=%.0f, \\zeta=%.0f, \\alpha=%.2f'...
             '\nGrid=%dx%d, R_\\infty=%.0e'], ...
    sol.Du, sol.zeta, sol.alpha, numel(g.r), numel(g.t), g.r(end));
title(t)
xlabel('Electric field (\beta)'); 
ylabel('Steady-state velocity (V)')
legend('Numerical', 'Error from linear')
grid on
ylim([1e-4 1e1])

% Compute the non-linear component
subplot 122
s = betas(2)/betas(1);
I = 2:numel(betas);
b = betas(I);
e = V(I) - s*V(I-1);

% Apply linear regression to log-log graph.
J = 1:numel(b);
r = linreg(log10(b(J)), log10(e(J)));
loglog(b, e, '.', betas, (10^r.b)*betas.^r.a, 'r-.'); 
title('Nonlinear component')
xlabel('Electric field (\beta)'); 
ylabel('e[n] = V[n] - V[n-1]*\beta[n]/\beta[n-1]')
legend('Numerical', sprintf('Fit \\propto \\beta^{%.1f}', r.a))
grid on

newname = sprintf('[%dx%d]_%.0e_alpha=%.2f', numel(g.r), numel(g.t), g.r(end), sol.alpha);
fprintf('>> %s\n', newname)
print('-depsc2', [newname '.eps']);

s = load(matname);
save([newname '.mat'], '-struct', 's');

