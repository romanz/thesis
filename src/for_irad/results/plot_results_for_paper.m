function plot_results_for_paper

a = [0.0 0.3 0.5];
Du = [0.5 0.5 1];
z = [6 6 10];
R = [30 100 300];
U1 = [3.69314718055995, 3.69314718055995, 4.25752957407993];
U3 = [0.0263988658183029, 0.304452832842969, 0.966868914065528];

figure(1); clf; 
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',[0 0 1200 400]);

for k = 1:3
    subplot(1,3,k)
    v = [];
    for j = 2
        f = sprintf('alpha=%.1f_Du=%.1f_zeta=%.1f_Rmax=%.0f.mat', a(k), Du(k), z(k), R(j));
        s = load(f);
        V1 = U1(k) * s.b;
        V3 = U3(k) * s.b.^3;
        v = [v, abs(s.Vr1 - V1)];
    end
    loglog(s.b, v, '.', s.b, V3, '-');
    t = sprintf('\\alpha = %.1f, Du = %.1f, \\zeta = %.1f', a(k), Du(k), z(k));
    title(t)
    xlabel('electric field magnitude (\beta)')
    ylabel('cubic correction to linear steady-state velocity')
    legend('Numerical', 'Analytical', ...
        'Location', 'SouthEast')
    xlim([min(s.b) max(s.b)])
    ylim([1e-6 10])
    grid on;    
end
print -depsc2 results1

figure(2); clf; 
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',[0 0 1200 400]);

for k = 1:3
    subplot(1,3,k)
    v = [];
    for j = 2
        f = sprintf('for_irad/results/alpha=%.1f_Du=%.1f_zeta=%.1f_Rmax=%.0f.mat', a(k), Du(k), z(k), R(j));
        s = load(f);
        V1 = U1(k) * s.b;
        V3 = U3(k) * s.b.^3;
        v = [v, s.Vr1];
    end
    loglog(s.b, v, '.', s.b, V1 + V3, '-');
    t = sprintf('\\alpha = %.1f, Du = %.1f, \\zeta = %.1f', a(k), Du(k), z(k));
    title(t)
    xlabel('electric field magnitude (\beta)')
    ylabel('steady-state velocity')
    legend('Numerical', 'Analytical', ...
        'Location', 'SouthEast')
    xlim([min(s.b) max(s.b)])
    ylim([6e-2 100])
    grid on;    
end
print -depsc2 results2

