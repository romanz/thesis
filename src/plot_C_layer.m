function plot_C_layer
clear
load results_C(beta)_512x512_Rmax=10.mat
b = betas(1:end-1);
clf;
I = 1:numel(b);
R = [];
w = [];

clf
c = C(:, 1:2:end);
plot(r-1, c, '-', 'MarkerSize', 5, 'LineWidth', 2)    
xlim([0 1])
ylim([0 1])
xlabel('Radial distance from particle surface', 'FontSize', 16)
ylabel('Ion concentration (C)', 'FontSize', 16)
print -depsc2 boundary_layer_C1

clf
c = C(1, :);
logc = log(c);
I = 10:50;
r = linreg(b(I), logc(I));
semilogy(b, c, 'o', b(I), c(I), 'o', b, exp(b*r.a + r.b), '-', 'MarkerSize', 5, 'LineWidth', 2)
xlim([0 5])
ylim([.1 1])
xlabel('Electric field magnitude (\beta)', 'FontSize', 16)
ylabel('C(r=1, \theta=0)', 'FontSize', 16)
func = sprintf('\n%.3f e^{%.3f \\beta}', exp(r.b), r.a);
legend('Numerical results', 'Points to fit', ['Exponential fit ' func])
    
print -depsc2 boundary_layer_C2
