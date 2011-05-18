%% Grab all streamline plots
load
clf;
for k = 1:numel(sols)
    sol = sols{k};
    clf;
    hold on
    g = sol.grid.Psi;
    contour(g.X .* cos(g.Y), g.X .* sin(g.Y), sol.Psi, sol.Vinf * linspace(-1, 1));
    t = sol.grid.theta;
    fill(cos(t), sin(t), 'k');
    hold off
    axis equal;
    d = 5;
    xlim([-d d]);
    ylim([0 d])
    text(-4.5, 4.5, sprintf('Streamlines - \\beta = %f', betas(k)))
    drawnow;
    print('-dpng', sprintf('%02d', k))
end
