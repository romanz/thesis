function animate
    S = load('matlab');
    figure(1)
    set(gcf, 'PaperPositionMode', 'auto');
%     set(gcf, 'Position',[100 100 900 600]);
    V = attr(S.solutions, 'Vinf');
    Va = attr(S.solutions, 'Vinf_approx');
    betas = attr(S.solutions, 'beta');
    for k = 1:numel(S.solutions)
        sol = S.solutions{k};
        s = [['\beta = ' sprintf('%.3f', sol.beta)] ...
            sprintf('  ') ...
            ['V = '  sprintf('%.3f', sol.Vinf)]];
        clf; 
        subplot(2,2,[1 3])
        colormap(jet)
        loglog(betas, V, '-b', betas, Va, ':r', ...
            sol.beta, sol.Vinf, '.b', ...
            sol.beta, sol.Vinf_approx, '.r', ...
            sol.beta, brenner(sol), '.k', ...
            'MarkerSize', 20, 'LineWidth', 2.5)
        title('Electrokinetic Solver Results', 'FontSize', 16)
        xlabel('Electric field magnitude (\beta)', 'FontSize', 12)
        ylabel('Particle velocity magnitude (V)', 'FontSize', 12)
        legend({'Numerical Solver', 'Analytical Approximation'}, ...
            'Location', 'SouthEast', 'FontSize', 8, 'Box', 'on')
        axis([0.01 1 0.001 1])
        text(0.04, 0.4, s, 'FontSize', 12);
        subplot(2,2,2);
        hold on
        show(@contour, sol, 'Psi', 10, +1, sol.Vinf * linspace(-1, 1, 50));
        show(@contour, sol, 'Psi', 10, -1, sol.Vinf * linspace(-1, 1, 50));
        sphere();
        hold off; 
        axis equal; axis([-4 4 -3 3]); 
        title('Fluid streamlines', 'FontSize', 14)

        subplot(2,2,4);
        hold on
        for d = [-1 1]
            show(@surf, sol, 'C', 10, d, 'EdgeColor', 'none');
        end
        set(gca, 'CLim', [0 2])
        shading('interp')
        colorbar;
        view(2)
        sphere();
        hold off; 
        axis equal; axis([-4 4 -3 3]); 
        title('Salt Concentration', 'FontSize', 14)
        drawnow;
        print('-dpng', '-r95', sprintf('results/%03d', k-1))
    end
end

function sphere()
    t = linspace(0, 2*pi, 100);
    fill(cos(t), sin(t), 'k');
end
