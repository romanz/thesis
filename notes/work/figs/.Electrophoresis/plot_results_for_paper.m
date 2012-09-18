function plot_results_for_paper

a = [0.0 0.3 0.5];
Du = [0.5 0.5 1];
z = [6 6 10];
R = [30 100 300];
U1 = [3.69314718055995, 3.69314718055995, 4.25752957407993];
U3 = [0.0263988658183029, 0.304452832842969, 0.966868914065528];

close all;

for k = 1:3
    for j = 2
        f = sprintf('alpha=%.1f_Du=%.1f_zeta=%.1f_Rmax=%.0f', a(k), Du(k), z(k), R(j));
        s = load([f '.mat']);

        b = linspace(min(s.b), max(s.b), 1e3);
        u = s.Vr1; % numerical solution
        e = abs(s.Vr1 - U1(k) * s.b); % deviation from linear regime

        figure(); hold on;        
        plot(b, U3(k) * b.^3, '-', 'LineWidth', 2, 'Color', 'k')
        plot(s.b, e,          'o', 'MarkerSize', 8, 'Color', 'k');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        t = sprintf('\\alpha = %.1f, Du = %.1f, \\zeta = %.1f', a(k), Du(k), z(k));
%         title(t, 'FontSize', 24)
        set(gca, 'FontSize', 16)
        xlabel('Electric field magnitude (\beta)', 'FontSize', 20)
        ylabel('Departure from linear analytic solution', 'FontSize', 20)
        xlim([3e-2 4])
        ylim([3e-5 10])
        legend({'Analytical', 'Numerical'}, 'Location', 'SouthEast', 'FontSize', 18);
        legend boxoff;
        grid on;    
        print('-depsc2', sprintf('departure_from_linear_%s.eps', f))

        figure(); hold on;        
        plot(b, U1(k) * b + U3(k) * b.^3, '-', 'LineWidth', 2, 'Color', 'k');
        plot(s.b, u,                      'o', 'MarkerSize', 8, 'Color', 'k');
        set(gca, 'XScale', 'log', 'YScale', 'log');
        t = sprintf('\\alpha = %.1f, Du = %.1f, \\zeta = %.1f', a(k), Du(k), z(k));
%         title(t, 'FontSize', 24)
        set(gca, 'FontSize', 16)
        xlabel('Electric field magnitude (\beta)', 'FontSize', 20)
        ylabel('Steady-state velocity', 'FontSize', 20)
        xlim([2e-2 4])
        ylim([8e-2 50])
        legend({'Analytical', 'Numerical'}, 'Location', 'SouthEast', 'FontSize', 18);
        legend boxoff;
        grid on;    
        print('-depsc2', sprintf('numerical_vs_analytical_%s.eps', f))
    end
end
