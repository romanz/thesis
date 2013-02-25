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

hold on;
I = 1:2:numel(q);
% for i = I
    f = @(x) x;
    plot(f(r-1), F(:, I), '-', 'LineWidth', 2);
    plot(f(q(I)-1), 0.5, 'o', 'MarkerSize', 5);
    xlim([0 0.2])
    ylim([0 1])
    xlabel('Radial distance from particle surface', 'FontSize', 16)
    ylabel('Normalized electric field radial component', 'FontSize', 16)
% end
print -depsc2 BoundaryLayerWidth_E1.eps

clf;
hold on;
I = 1:2:numel(q);
% for i = I
    plot((r-1) * exp(b(I) / 2), F(:, I), '-', 'LineWidth', 2);
    plot((q(I)-1) .* exp(b(I) / 2), 0.5, 'o', 'MarkerSize', 5);
    xlim([0 1])
    ylim([0 1])
    xlabel('Radial distance from particle surface (scaled by e^{\beta/2})', 'FontSize', 16)
    ylabel('Normalized electric field radial component', 'FontSize', 16)
% end
print -depsc2 BoundaryLayerWidth_E2.eps

clf;
%hold on
y = log(q-1);
I = 10:50;
r = linreg(b(I), y(I))
semilogy(b, q-1, 'o', b(I), q(I)-1, 'o', b, exp(r.a*b + r.b), '-' ...
    , 'MarkerSize', 5, 'LineWidth', 2)
xlim([0 5])
ylim([0.01 0.2])
func = sprintf('\n%.3f e^{%.3f \\beta}', exp(r.b), r.a);
legend({'Numerical results', 'Points to fit', ['Exponential fit ' func]}, ...
    'FontSize', 16)
xlabel('\beta', 'FontSize', 16)
ylabel('E_r boudary layer width', 'FontSize', 16)
print -depsc2 BoundaryLayerWidth_E3.eps

clf; 
e_ = E(1, :) ./ b;
I = 10:50;
r = linreg(b(I), log(e_(I)))
semilogy(b, e_, 'o', b(I), e_(I), 'o', b, exp(r.a * b + r.b), '-' ...
        , 'MarkerSize', 5, 'LineWidth', 2)
xlabel('\beta', 'FontSize', 16)
ylabel('E_r(r=1, \theta=0) / \beta', 'FontSize', 16)
xlim([0 5])
ylim([1 10])
func = sprintf('\n%.3f e^{%.3f \\beta}', exp(r.b), r.a);
legend('Numerical results', 'Points to fit', ['Exponential fit ' func], ...
    'Location', 'SouthEast')
print -depsc2 BoundaryLayerWidth_E4.eps
end