function plot_C_Er_product

load results_Er(beta)_512x512_Rmax=10.mat
load results_C(beta)_512x512_Rmax=10.mat

b = betas(1:end-1);

plot(b, E(1, :) .* C(1, :) ./ b, 'LineWidth', 2)
xlim([0 5])
xlabel('\beta', 'FontSize', 16)
ylabel('C \cdot E_r / \beta at r=1, \theta=0', 'FontSize', 16, 'FontSize', 16)
print -depsc2 C_Er_product