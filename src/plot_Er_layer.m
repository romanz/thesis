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
subplot(1,2,1)
hold on;
I = 1:2:numel(q);
plot(r-1, F(:, I), '-');
plot(q(I)-1, 0.5, '.');
xlim([0 0.2])
ylim([0 1])
xlabel('Radial distance from particle surface')
ylabel('Normalized electric field radial component')

subplot(1,2,2)
%hold on
y = log(q-1);
I = 2:2:68;
r = linreg(b(I), y(I))
semilogy(b(I), q(I)-1, '.', b, exp(r.a*b + r.b))
ylim([0.005 0.2])
legend('Numerical results', 'Exponential fit')
xlabel('\beta')
ylabel('Boudary layer width')

print -depsc2 BoundaryLayerWidth.eps
end