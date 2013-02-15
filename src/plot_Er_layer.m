function plot_Er_layer
clear
load results_Er(beta)_512x512_Rmax=10.mat
b = betas(1:end-1);
F = (E - repmat(E(end, :), size(E, 1), 1));
F = F * diag(1./max(F));
for j = 1:size(F, 2)
    f = F(:, j);
    q(j) = interp1(f, r, 0.5);
end
clf;
subplot(2,2,1)
hold on;
I = 1:2:numel(q);
for i = I
    f = @(x) x;
    plot(f(r-1), F(:, i), '-');
    plot(f(q(i)-1), 0.5, '.');
    drawnow;
    pause(.1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Radial distance from particle surface')
    ylabel('Normalized electric field radial component')
end
subplot(2,2,2)
hold on;
I = 1:2:numel(q);
for i = I
    f = @(x) x * exp(b(i) / 2);
    plot(f(r-1), F(:, i), '-');
    plot(f(q(i)-1), 0.5, '.');
    drawnow;
    pause(.1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Radial distance from particle surface (scaled by e^{\beta/2})')
    ylabel('Normalized electric field radial component')
end
subplot(2,2,3)
%hold on
y = log(q-1);
I = 10:50;
r = linreg(b(I), y(I))
semilogy(b, q-1, '.', b(I), q(I)-1, '.', b, exp(r.a*b + r.b))
ylim([0.005 0.2])
func = sprintf('\n%.3f e^{%.3f \\beta}', exp(r.b), r.a);
legend('Numerical results', 'Points to fit', ['Exponential fit ' func])
xlabel('\beta')
ylabel('Boudary layer width')

subplot(2,2,4)
%hold on
semilogy(b, E(1, :) ./ b, '.')
xlabel('\beta')
ylabel('E_r(r=1, \theta=0)')

print -depsc2 BoundaryLayerWidth.eps
end