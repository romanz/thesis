function plot_C_layer
clear
load results_C(beta)_512x512_Rmax=10.mat
b = betas(1:end-1);
clf;
I = 1:numel(b);
R = [];
w = [];
for i = I
    subplot(2,2,1)
    hold on;
    c = C(:, i);
    plot(r-1, c)    
    xlim([0 1])
    ylim([0 1])
    xlabel('Radial distance from particle surface')
    ylabel('Ion concentration (C)')
    
    subplot(2,2,2)
    hold on;
    c_ = c - c(1);
    c_ = c_ / c_(end);
    I = c_ < 0.999;
    w(i) = interp1(c_(I), r(I)-1, 0.5);
    plot(r-1, c_, w, 0.5, '.k')
    xlim([0 1])
    ylim([0 1])
    xlabel('Radial distance from particle surface')
    ylabel('Normalized Ion concentration')
    
    subplot(2,2,3)
    semilogy(b(1:i), C(1, 1:i), '.-')
    xlim([0 7])
    ylim([3e-2 1])
    xlabel('Electric field magnitude (\beta)')
    ylabel('C(r=1, \theta=0)')
    
    subplot(2,2,4)
    loglog(b(1:i), w(1:i), '.-')
    xlim([0.1 7])
    ylim([0.08 0.5])
    xlabel('Electric field magnitude (\beta)')
    ylabel('Boundary layer width')
    
    drawnow;
    pause(.1)
end
return
subplot(1,2,2)
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

print -depsc2 BoundaryLayerWidth.eps
% end