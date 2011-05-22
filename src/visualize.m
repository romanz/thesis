function visualize
    load
    clf;
    for k = 1:numel(sols)
        sol = sols{k};
        clf; 
        hold on
        show(@contour, sol, 'Psi', 10, sol.Vinf * linspace(-1, 1));
        sphere(sol);
%         show(@surf, sol, 'Phi', 3, 'EdgeColor', 'k');
        hold off; axis equal; axis([-1 1 0 1]*5)
        
        title(sprintf('\\beta = %f : V_{\\infty} = %f', sol.beta, sol.Vinf))
        drawnow;
        print('-dpng', sprintf('%02d', k))
    end
end

function sphere(sol)
    t = sol.grid.theta;
    fill(cos(t), sin(t), 'k');
end